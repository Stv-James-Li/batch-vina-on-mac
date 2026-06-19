#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import numpy as np

def parse_pdbqt_coordinates(pdbqt_file):
    """
    读取PDBQT文件，提取所有原子坐标，并计算中心和大小
    
    Args:
        pdbqt_file: PDBQT文件路径
        
    Returns:
        center_x, center_y, center_z: 分子的中心坐标
        size_x, size_y, size_z: 分子在三个维度上的大小
    """
    coordinates = []
    
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    # PDBQT格式中原子坐标通常在第7-9列
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates.append([x, y, z])
                except ValueError:
                    continue
    
    if not coordinates:
        raise ValueError(f"No coordinates found in {pdbqt_file}")
    
    coordinates = np.array(coordinates)
    
    # 计算分子的边界盒子
    min_coords = np.min(coordinates, axis=0)
    max_coords = np.max(coordinates, axis=0)
    
    # 计算中心坐标
    center_x = (min_coords[0] + max_coords[0]) / 2
    center_y = (min_coords[1] + max_coords[1]) / 2
    center_z = (min_coords[2] + max_coords[2]) / 2
    
    # 计算大小（添加一些余量可能会更好）
    size_x = max_coords[0] - min_coords[0]
    size_y = max_coords[1] - min_coords[1]
    size_z = max_coords[2] - min_coords[2]
    
    return center_x, center_y, center_z, size_x, size_y, size_z

def generate_vina_config(pdbqt_file, output_dir=None):
    """
    为PDBQT文件生成Vina配置文件
    
    Args:
        pdbqt_file: PDBQT文件路径
        output_dir: 输出目录，默认为PDBQT文件所在目录
    """
    try:
        # 提取分子中心和大小
        center_x, center_y, center_z, size_x, size_y, size_z = parse_pdbqt_coordinates(pdbqt_file)
        
        # 准备输出文件名
        base_name = os.path.splitext(os.path.basename(pdbqt_file))[0]
        if output_dir is None:
            output_dir = os.path.dirname(pdbqt_file)
        output_file = os.path.join(output_dir, f"{base_name}.txt")
        
        # 写入配置文件
        with open(output_file, 'w') as f:
            f.write(f"center_x={center_x:.3f}\n")
            f.write(f"center_y={center_y:.3f}\n")
            f.write(f"center_z={center_z:.3f}\n")
            f.write(f"size_x={size_x:.3f}\n")
            f.write(f"size_y={size_y:.3f}\n")
            f.write(f"size_z={size_z:.3f}\n")
            f.write("num_modes=9\n")
        
        print(f"生成配置文件: {output_file}")
        return True
    except Exception as e:
        print(f"处理 {pdbqt_file} 时出错: {str(e)}")
        return False

def process_directory(directory_path, output_dir=None):
    """
    处理指定目录下的所有PDBQT文件
    
    Args:
        directory_path: 包含PDBQT文件的目录路径
        output_dir: 配置文件输出目录，默认与PDBQT文件相同
    """
    # 确保目录路径有效
    if not os.path.isdir(directory_path):
        print(f"错误: {directory_path} 不是有效目录")
        return
    
    # 如果指定了输出目录，确保它存在
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 查找所有PDBQT文件
    pdbqt_files = glob.glob(os.path.join(directory_path, "*.pdbqt"))
    
    if not pdbqt_files:
        print(f"在 {directory_path} 中未找到PDBQT文件")
        return
    
    success_count = 0
    for pdbqt_file in pdbqt_files:
        if generate_vina_config(pdbqt_file, output_dir):
            success_count += 1
    
    print(f"处理完成。成功生成 {success_count}/{len(pdbqt_files)} 个配置文件。")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='从PDBQT文件生成Vina配置文件')
    parser.add_argument('--directory', '-d', default='.', help='包含PDBQT文件的目录路径（默认为当前目录）')
    parser.add_argument('--output', '-o', help='配置文件输出目录（默认与PDBQT文件相同）')
    
    args = parser.parse_args()
    
    process_directory(args.directory, args.output)