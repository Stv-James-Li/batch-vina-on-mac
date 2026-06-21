#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

PLIP_STYLES = {
    'hydrogen_bond': ('HBonds', 'yellow', '2.4'),
    'halogen_bond': ('HalogenBonds', 'magenta', '2.6'),
    'hydrophobic_interaction': ('Hydrophobic', 'gray50', '1.8'),
    'salt_bridge': ('SaltBridges', 'tv_blue', '2.6'),
    'pi_stack': ('PiStacks', 'orange', '2.4'),
    'pi_cation_interaction': ('PiCation', 'purple', '2.4'),
    'water_bridge': ('WaterBridges', 'cyan', '2.0'),
    'metal_complex': ('MetalComplexes', 'slate', '2.4'),
}

PLIP_GROUP_TO_ITEM = {
    'hydrogen_bonds': 'hydrogen_bond',
    'halogen_bonds': 'halogen_bond',
    'hydrophobic_interactions': 'hydrophobic_interaction',
    'salt_bridges': 'salt_bridge',
    'pi_stacks': 'pi_stack',
    'pi_cation_interactions': 'pi_cation_interaction',
    'water_bridges': 'water_bridge',
    'metal_complexes': 'metal_complex',
}


def run_in_autodock(command: str):
    subprocess.run([
        'bash', '-lc',
        'source /Users/jamesli/miniforge3/etc/profile.d/conda.sh && conda activate autodock >/dev/null && ' + command,
    ], check=True)


def quote_path(path: Path) -> str:
    return str(path).replace('\\', '\\\\').replace('"', '\\"')


def first_model_lines(path: Path):
    out = []
    seen = False
    in_model = False
    for line in path.read_text(errors='ignore').splitlines():
        if line.startswith('MODEL'):
            if seen:
                break
            seen = True
            in_model = True
            continue
        if line.startswith('ENDMDL') and in_model:
            break
        if in_model or not seen:
            out.append(line)
    return out


def element_from_line(line: str) -> str:
    elem = line[76:78].strip() if len(line) >= 78 else ''
    if elem:
        return elem.upper()
    tail = line[66:].strip().split()
    if tail:
        atom_type = ''.join(ch for ch in tail[-1] if ch.isalpha()).upper()
        atom_type = {'OA': 'O', 'HD': 'H', 'A': 'C', 'NA': 'N', 'SA': 'S'}.get(atom_type, atom_type)
        if atom_type:
            return atom_type[:2] if atom_type[:2] in ('CL', 'BR') else atom_type[:1]
    name = line[12:16].strip().upper()
    return name[:2] if name[:2] in ('CL', 'BR') else name[:1]


def ligand_atom_name(element: str, counts: dict[str, int]) -> str:
    element = element.upper()
    counts[element] = counts.get(element, 0) + 1
    return f'{element}{counts[element]}'[:4]


def parse_xyz(line: str):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def format_pdb_line(record, serial, atom_name, resname, chain, resseq, x, y, z, element, occ=1.0, bfac=0.0):
    atom_name = atom_name[:4]
    atom_field = f' {atom_name:<3}' if len(atom_name) < 4 and atom_name[:1].isalpha() else f'{atom_name:<4}'
    return (
        f'{record:<6}{serial:>5} {atom_field} {resname[:3]:>3} {chain[:1]}{int(resseq):>4}    '
        f'{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{bfac:>6.2f}          {element[:2]:>2}\n'
    )


def receptor_to_pdb_lines(receptor_pdbqt: Path, chain='A'):
    lines = []
    serial = 1
    for line in first_model_lines(receptor_pdbqt):
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        resname = line[17:20].strip() or 'UNK'
        try:
            resseq = int(line[22:26])
        except ValueError:
            resseq = 1
        x, y, z = parse_xyz(line)
        element = element_from_line(line)
        atom_name = line[12:16].strip() or element
        lines.append(format_pdb_line('ATOM', serial, atom_name, resname, chain, resseq, x, y, z, element))
        serial += 1
    return lines, serial


def obabel_ligand_to_pdb(ligand_pdbqt: Path, output_pdb: Path):
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    run_in_autodock(f'obabel -ipdbqt "{quote_path(ligand_pdbqt)}" -opdb -O "{quote_path(output_pdb)}" -f 1 -l 1')


def read_obabel_ligand(obabel_pdb: Path, start_serial: int, chain='B', resname='LIG', resseq=1):
    atom_lines = []
    conect_lines = []
    old_to_new = {}
    serial = start_serial
    name_counts = {}
    source_lines = obabel_pdb.read_text(errors='ignore').splitlines()

    for line in source_lines:
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        old_serial = int(line[6:11])
        x, y, z = parse_xyz(line)
        element = element_from_line(line)
        atom_name = ligand_atom_name(element, name_counts)
        old_to_new[old_serial] = serial
        atom_lines.append(format_pdb_line('HETATM', serial, atom_name, resname, chain, resseq, x, y, z, element))
        serial += 1

    for line in source_lines:
        if not line.startswith('CONECT'):
            continue
        nums = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
        if not nums or nums[0] not in old_to_new:
            continue
        mapped = [old_to_new[n] for n in nums if n in old_to_new]
        if len(mapped) > 1:
            conect_lines.append('CONECT' + ''.join(f'{n:>5}' for n in mapped) + '\n')
    return atom_lines, conect_lines


def write_complex(receptor_lines, ligand_lines, conect_lines, output_pdb: Path):
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text(''.join(receptor_lines) + 'TER\n' + ''.join(ligand_lines) + ''.join(conect_lines) + 'TER\nEND\n')


def run_plip(complex_pdb: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    run_in_autodock(f'plip -f "{quote_path(complex_pdb)}" -o "{quote_path(output_dir)}" -x -t --pymol')


def xml_text(node, tag, default=''):
    value = node.findtext(tag)
    return value.strip() if value else default


def xml_coord(node, tag):
    elem = node.find(tag)
    if elem is None:
        return None
    try:
        return tuple(float(elem.findtext(axis)) for axis in ('x', 'y', 'z'))
    except (TypeError, ValueError):
        return None


def parse_plip_interactions(xml_path: Path):
    root = ET.parse(xml_path).getroot()
    interactions = []
    for group in root.findall('.//interactions/*'):
        item_tag = PLIP_GROUP_TO_ITEM.get(group.tag)
        if not item_tag:
            continue
        for item in group.findall(item_tag):
            lig = xml_coord(item, 'ligcoo')
            prot = xml_coord(item, 'protcoo')
            if lig is None or prot is None:
                continue
            interactions.append({
                'type': item_tag,
                'lig': lig,
                'prot': prot,
                'resnr': xml_text(item, 'resnr'),
                'restype': xml_text(item, 'restype'),
                'reschain': xml_text(item, 'reschain'),
            })
    return interactions


def write_plip_style_pml(structure_pdb: Path, xml_path: Path, output_pse: Path, pml_path: Path):
    interactions = parse_plip_interactions(xml_path)
    lines = [
        'reinitialize\n',
        f'load "{quote_path(structure_pdb)}", complex\n',
        'remove solvent\n',
        'select lig_sel, complex and resn LIG and chain B\n',
        'select prot_sel, complex and polymer.protein\n',
        'create protein, prot_sel\n',
        'create ligand, lig_sel\n',
        'disable complex\n',
        'hide everything, all\n',
        'bg_color white\n',
        'set opaque_background, on\n',
        'show cartoon, protein\n',
        'color gray80, protein\n',
        'show sticks, ligand\n',
        'show spheres, ligand and elem F+O+S+H\n',
        'set sphere_scale, 0.22, ligand\n',
        'color cyan, ligand and elem C\n',
        'color red, ligand and elem O\n',
        'color yellow, ligand and elem S\n',
        'color green, ligand and elem F\n',
        'color white, ligand and elem H\n',
        'show sticks, byres (protein within 4 of ligand)\n',
        'color marine, byres (protein within 4 of ligand)\n',
        'set stick_radius, 0.18\n',
        'set valence, 1\n',
        'set dash_gap, 0.25\n',
        'set dash_radius, 0.08\n',
        'set label_color, black\n',
        'set label_size, 18\n',
        'set label_outline_color, white\n',
    ]

    contact_residues = {}
    for interaction in interactions:
        contact_residues[(interaction['reschain'], interaction['resnr'])] = interaction['restype']
    for (chain, resnr), restype in sorted(contact_residues.items(), key=lambda x: (x[0][0], int(x[0][1]) if x[0][1].isdigit() else x[0][1])):
        label_name = f'label_{restype}{resnr}_{chain}'
        lines.extend([
            f'select {label_name}, protein and chain {chain} and resi {resnr} and name CA\n',
            f'label {label_name}, "{restype}{resnr}"\n',
            f'show sticks, protein and chain {chain} and resi {resnr}\n',
            f'color marine, protein and chain {chain} and resi {resnr}\n',
        ])

    counts = {}
    for interaction in interactions:
        group, color, width = PLIP_STYLES[interaction['type']]
        counts[group] = counts.get(group, 0) + 1
        name = f"{group}_{counts[group]}_{interaction['restype']}{interaction['resnr']}"
        lig_name = f'pt_{name}_lig'
        prot_name = f'pt_{name}_prot'
        lx, ly, lz = interaction['lig']
        px, py, pz = interaction['prot']
        lines.extend([
            f'pseudoatom {lig_name}, pos=[{lx:.3f},{ly:.3f},{lz:.3f}]\n',
            f'pseudoatom {prot_name}, pos=[{px:.3f},{py:.3f},{pz:.3f}]\n',
            f'distance {name}, {lig_name}, {prot_name}\n',
            f'color {color}, {name}\n',
            f'set dash_width, {width}, {name}\n',
            f'hide labels, {name}\n',
            f'hide everything, {lig_name} or {prot_name}\n',
        ])
    for group, _, _ in set(PLIP_STYLES.values()):
        lines.append(f'group {group}, {group}_*\n')
    lines.extend(['orient ligand\n', 'zoom ligand, 8\n', f'cmd.save("{quote_path(output_pse)}")\n'])
    pml_path.write_text(''.join(lines))
    return counts


def generate_plip_style_pse(plip_dir: Path, output_dir: Path, ligand_base: str):
    xml_path = plip_dir / f'{ligand_base}_complex_report.xml'
    structure_pdb = plip_dir / f'{ligand_base}_complex_protonated.pdb'
    if not xml_path.exists() or not structure_pdb.exists():
        print(f'  ! PLIP-style PSE skipped; missing XML/protonated PDB in {plip_dir}')
        return False
    output_dir.mkdir(parents=True, exist_ok=True)
    output_pse = output_dir / f'{ligand_base}_plip_style_labeled.pse'
    pml_path = output_pse.with_suffix('.pml')
    counts = write_plip_style_pml(structure_pdb, xml_path, output_pse, pml_path)
    run_in_autodock(f'pymol -cq "{quote_path(pml_path)}"')
    summary = ', '.join(f'{name}={count}' for name, count in sorted(counts.items())) or 'no interactions'
    print(f'  ✓ PLIP-style PSE: {output_pse.name} ({summary})')
    return True


def batch_create_complexes(receptor_folder, results_folder, output_folder='complex_pdbs_plip', run_plip_analysis=True, plip_folder=None, make_style_pse=True):
    receptor_folder = Path(receptor_folder)
    results_folder = Path(results_folder)
    output_folder = Path(output_folder)
    ligand_cache = output_folder / 'ligands_obabel'
    plip_folder = output_folder / 'plip' if plip_folder is None else Path(plip_folder)
    style_folder = output_folder / 'plip_style_pse_labeled'

    receptor_files = sorted(p for p in receptor_folder.glob('*.pdbqt') if not p.name.endswith('.pdbqt.pdbqt'))
    for receptor_path in receptor_files:
        receptor_base = receptor_path.stem
        receptor_lines, ligand_start = receptor_to_pdb_lines(receptor_path)
        ligands = sorted(results_folder.glob(f'{receptor_base}_*.pdbqt'))
        if not ligands:
            print(f'[SKIP] No ligands for {receptor_base}')
            continue
        print(f'[{receptor_base}] → {len(ligands)} ligand(s)')
        for ligand_path in ligands:
            ligand_base = ligand_path.stem
            obabel_pdb = ligand_cache / f'{ligand_base}_obabel.pdb'
            complex_pdb = output_folder / f'{ligand_base}_complex.pdb'
            obabel_ligand_to_pdb(ligand_path, obabel_pdb)
            ligand_lines, conect_lines = read_obabel_ligand(obabel_pdb, ligand_start)
            write_complex(receptor_lines, ligand_lines, conect_lines, complex_pdb)
            print(f'  ✓ {complex_pdb.name} (OpenBabel ligand bonds)')
            if run_plip_analysis:
                current_plip_dir = plip_folder / ligand_base
                run_plip(complex_pdb, current_plip_dir)
                if make_style_pse:
                    generate_plip_style_pse(current_plip_dir, style_folder, ligand_base)
    print('\nAll done.')


def main():
    parser = argparse.ArgumentParser(description='用 OpenBabel 组合受体-配体复合物，运行 PLIP，并生成白底带接触残基标签的 PLIP-style PSE。')
    parser.add_argument('--receptor-folder', default='receptor')
    parser.add_argument('--results-folder', default='results')
    parser.add_argument('--output-folder', default='complex_pdbs_plip')
    parser.add_argument('--plip-folder', default=None)
    parser.add_argument('--no-plip', action='store_true')
    parser.add_argument('--no-style-pse', action='store_true', help='不生成白底、带残基标签的 PLIP-style PSE')
    args = parser.parse_args()
    batch_create_complexes(
        receptor_folder=args.receptor_folder,
        results_folder=args.results_folder,
        output_folder=args.output_folder,
        run_plip_analysis=not args.no_plip,
        plip_folder=args.plip_folder,
        make_style_pse=not args.no_style_pse,
    )


if __name__ == '__main__':
    main()
