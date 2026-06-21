import argparse
import os
import shutil
import subprocess


def pdbqt_to_pdb_lines(pdbqt_path):
    """
    从PDBQT文件中提取PDB格式的原子行，保留原始原子序号和链ID。
    仅保留ATOM/HETATM/TER/CONECT行，去掉PDBQT特有的ROOT/BRANCH/TORSDOF等行。
    """
    pdb_lines = []
    with open(pdbqt_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record in ("ATOM", "HETATM", "TER", "CONECT", "REMARK", "MODEL", "ENDMDL"):
                pdb_lines.append(line[:80].rstrip() + "\n")
    return pdb_lines


def extract_first_model(lines):
    """如果PDBQT含多个MODEL（多构象），只取第一个MODEL的原子行。"""
    result = []
    in_model = False
    model_found = False
    for line in lines:
        record = line[:6].strip()
        if record == "MODEL":
            if model_found:
                break
            in_model = True
            model_found = True
            continue
        if record == "ENDMDL" and in_model:
            break
        if in_model or not model_found:
            result.append(line)
    return result


def format_pdb_atom_line(
    record,
    serial,
    atom_name,
    res_name,
    chain_id,
    res_seq,
    x,
    y,
    z,
    occupancy=1.0,
    bfactor=0.0,
    element="",
):
    """按标准PDB定宽列写出ATOM/HETATM行。"""
    atom_name = atom_name[:4]
    if len(atom_name.strip()) < 4 and atom_name[:1].isalpha():
        atom_field = f" {atom_name.strip():<3}"
    else:
        atom_field = f"{atom_name:<4}"
    return (
        f"{record:<6}{serial:>5} "
        f"{atom_field} "
        f"{res_name[:3]:>3} {chain_id[:1]}{int(res_seq):>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{occupancy:>6.2f}{bfactor:>6.2f}          "
        f"{element[:2]:>2}\n"
    )


def parse_pdbqt_atom_fields(line):
    """从PDBQT/PDB原子行提取写标准PDB需要的字段。"""
    atom_name = line[12:16]
    res_name = line[17:20].strip() or "UNK"
    chain_id = (line[21].strip() or "A") if len(line) > 21 else "A"
    try:
        res_seq = int(line[22:26])
    except ValueError:
        res_seq = 1
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    try:
        occupancy = float(line[54:60])
    except ValueError:
        occupancy = 1.0
    try:
        bfactor = float(line[60:66])
    except ValueError:
        bfactor = 0.0
    element = line[76:78].strip() if len(line) >= 78 else ""
    if not element:
        pdbqt_type = line[77:].strip().split()[-1] if line[77:].strip().split() else ""
        element = "".join(ch for ch in pdbqt_type if ch.isalpha())[:2]
    if not element:
        element = atom_name.strip()[:1]
    return atom_name, res_name, chain_id, res_seq, x, y, z, occupancy, bfactor, element


def renumber_atoms(atom_lines, start_serial=1, chain_id=None, hetatm_residue_name=None, residue_number=None):
    """
    重新连续编号原子序号（第7-11列），可选：强制指定链ID、残基名。
    返回 (重编号后的行列表, 下一个可用的serial编号, 旧serial→新serial的映射)
    """
    new_lines = []
    serial = start_serial
    serial_map = {}

    for line in atom_lines:
        record = line[:6].strip()
        if record in ("ATOM", "HETATM"):
            old_serial_str = line[6:11].strip()
            try:
                old_serial = int(old_serial_str)
            except ValueError:
                old_serial = None

            atom_name, res_name, original_chain, res_seq, x, y, z, occupancy, bfactor, element = parse_pdbqt_atom_fields(line)
            if hetatm_residue_name and record == "HETATM":
                res_name = hetatm_residue_name
            if chain_id:
                original_chain = chain_id
            if residue_number is not None:
                res_seq = residue_number

            new_lines.append(format_pdb_atom_line(
                record=record,
                serial=serial,
                atom_name=atom_name,
                res_name=res_name,
                chain_id=original_chain,
                res_seq=res_seq,
                x=x,
                y=y,
                z=z,
                occupancy=occupancy,
                bfactor=bfactor,
                element=element,
            ))

            if old_serial is not None:
                serial_map[old_serial] = serial
            serial += 1
        elif record == "TER":
            new_lines.append(f"TER   {serial:>5}\n")
            serial += 1
        elif record == "CONECT":
            new_lines.append(("CONECT_PENDING", line))
        else:
            new_lines.append(line)

    return new_lines, serial, serial_map


def remap_conect_lines(pending_lines, serial_map):
    """根据serial_map将CONECT行中的旧原子序号替换为新序号。"""
    result = []
    for item in pending_lines:
        if isinstance(item, tuple) and item[0] == "CONECT_PENDING":
            line = item[1]
            nums = []
            for i in range(1, 5):
                start = 6 + (i - 1) * 5
                end = start + 5
                seg = line[start:end].strip()
                if seg:
                    try:
                        old_n = int(seg)
                        nums.append(serial_map.get(old_n, old_n))
                    except ValueError:
                        pass
            if nums:
                result.append(f"CONECT{''.join(f'{n:>5}' for n in nums)}\n")
        else:
            result.append(item)
    return result


def merge_receptor_ligand_to_pdb(
    receptor_pdbqt,
    ligand_pdbqt,
    output_pdb,
    receptor_chain="A",
    ligand_chain="B",
    ligand_resname="LIG",
    ligand_resnum=1,
):
    """将受体和配体PDBQT合并为规范的复合物PDB文件。"""
    rec_raw = extract_first_model(pdbqt_to_pdb_lines(receptor_pdbqt))
    rec_atom_lines = [l for l in rec_raw if l[:6].strip() in ("ATOM", "HETATM", "TER")]
    rec_numbered, next_serial, rec_serial_map = renumber_atoms(
        rec_atom_lines,
        start_serial=1,
        chain_id=receptor_chain,
    )
    rec_numbered = remap_conect_lines(rec_numbered, rec_serial_map)

    lig_raw = extract_first_model(pdbqt_to_pdb_lines(ligand_pdbqt))
    lig_atom_lines = []
    for line in lig_raw:
        record = line[:6].strip()
        if record in ("ATOM", "HETATM"):
            lig_atom_lines.append("HETATM" + line[6:])
        elif record == "CONECT":
            lig_atom_lines.append(("CONECT_PENDING", line))

    lig_numbered, _, lig_serial_map = renumber_atoms(
        lig_atom_lines,
        start_serial=next_serial,
        chain_id=ligand_chain,
        hetatm_residue_name=ligand_resname,
        residue_number=ligand_resnum,
    )
    lig_numbered = remap_conect_lines(lig_numbered, lig_serial_map)

    os.makedirs(os.path.dirname(output_pdb) or ".", exist_ok=True)
    with open(output_pdb, "w") as out:
        out.write(f"REMARK  Complex: {os.path.basename(receptor_pdbqt)} + {os.path.basename(ligand_pdbqt)}\n")
        out.write(f"REMARK  Receptor chain: {receptor_chain} | Ligand chain: {ligand_chain} | Ligand resname: {ligand_resname}\n")
        for line in rec_numbered:
            if isinstance(line, str):
                out.write(line)
        out.write("TER\n")
        for line in lig_numbered:
            if isinstance(line, str):
                out.write(line)
        out.write("TER\nEND\n")

    print(f"  ✓ {os.path.basename(output_pdb)}")


def run_plip_analysis(complex_pdb, plip_root):
    """对单个复合物PDB运行PLIP，并尽量输出XML/TXT报告和PyMOL .pse会话。"""
    if shutil.which("plip") is None:
        print("  ! 未找到 plip 命令，跳过PLIP分析。安装：pip install plip")
        return False

    complex_name = os.path.splitext(os.path.basename(complex_pdb))[0]
    output_dir = os.path.join(plip_root, complex_name)
    os.makedirs(output_dir, exist_ok=True)

    commands = [
        ["plip", "-f", complex_pdb, "-o", output_dir, "-x", "-t", "-y"],
        ["plip", "-f", complex_pdb, "-o", output_dir, "-x", "-t", "--pymol"],
    ]

    last_error = None
    for command in commands:
        try:
            subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            pse_files = [f for f in os.listdir(output_dir) if f.endswith(".pse")]
            if pse_files:
                print(f"  ✓ PLIP: {output_dir} ({len(pse_files)} pse)")
                return True

            last_error = "PLIP运行成功但未生成.pse文件"
        except subprocess.CalledProcessError as exc:
            last_error = exc.stderr.strip() or exc.stdout.strip() or str(exc)

    report_command = ["plip", "-f", complex_pdb, "-o", output_dir, "-x", "-t", "-p"]
    try:
        subprocess.run(report_command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"  ! PLIP已输出报告/图片，但未生成pse: {output_dir}。原因：{last_error}")
        return False
    except subprocess.CalledProcessError as exc:
        last_error = exc.stderr.strip() or exc.stdout.strip() or str(exc)

    print(f"  ✗ PLIP failed for {os.path.basename(complex_pdb)}: {last_error}")
    return False


def batch_create_complexes(
    receptor_folder,
    results_folder,
    output_folder="complex_pdbs_plip",
    receptor_chain="A",
    ligand_chain="B",
    ligand_resname="LIG",
    run_plip=True,
    plip_folder=None,
):
    os.makedirs(output_folder, exist_ok=True)
    if plip_folder is None:
        plip_folder = os.path.join(output_folder, "plip")

    receptor_files = sorted([
        f for f in os.listdir(receptor_folder)
        if f.endswith(".pdbqt") and not f.endswith(".pdbqt.pdbqt")
    ])

    for receptor_file in receptor_files:
        base_name = os.path.splitext(receptor_file)[0]
        legacy_base_name = f"{base_name}.pdbqt"
        receptor_path = os.path.join(receptor_folder, receptor_file)
        preferred_ligands = sorted([
            f for f in os.listdir(results_folder)
            if f.startswith(base_name + "_") and f.endswith(".pdbqt")
        ])
        legacy_ligands = sorted([
            f for f in os.listdir(results_folder)
            if f.startswith(legacy_base_name + "_") and f.endswith(".pdbqt")
        ])
        matching_ligands = preferred_ligands or legacy_ligands

        if not matching_ligands:
            print(f"[SKIP] No ligands for {base_name}")
            continue
        if legacy_ligands and not preferred_ligands:
            print(f"[WARN] Using legacy docking result names for {base_name}: {legacy_base_name}_*.pdbqt")

        print(f"[{base_name}] → {len(matching_ligands)} ligand(s)")
        for ligand_file in matching_ligands:
            ligand_base = os.path.splitext(ligand_file)[0]
            ligand_path = os.path.join(results_folder, ligand_file)
            output_pdb = os.path.join(output_folder, f"{ligand_base}_complex.pdb")

            try:
                merge_receptor_ligand_to_pdb(
                    receptor_pdbqt=receptor_path,
                    ligand_pdbqt=ligand_path,
                    output_pdb=output_pdb,
                    receptor_chain=receptor_chain,
                    ligand_chain=ligand_chain,
                    ligand_resname=ligand_resname,
                    ligand_resnum=1,
                )
                if run_plip:
                    run_plip_analysis(output_pdb, plip_folder)
            except Exception as e:
                print(f"  ✗ {ligand_file}: {e}")

    print("\nAll done.")


def main():
    parser = argparse.ArgumentParser(description="生成受体-配体复合物PDB，并对每个复合物运行PLIP分析。")
    parser.add_argument("--receptor-folder", default="receptor")
    parser.add_argument("--results-folder", default="results")
    parser.add_argument("--output-folder", default="complex_pdbs_plip")
    parser.add_argument("--plip-folder", default=None)
    parser.add_argument("--no-plip", action="store_true", help="只生成复合物PDB，不运行PLIP")
    args = parser.parse_args()

    batch_create_complexes(
        receptor_folder=args.receptor_folder,
        results_folder=args.results_folder,
        output_folder=args.output_folder,
        receptor_chain="A",
        ligand_chain="B",
        ligand_resname="LIG",
        run_plip=not args.no_plip,
        plip_folder=args.plip_folder,
    )


if __name__ == "__main__":
    main()
