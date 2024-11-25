#!/usr/bin/env python3
"""
Convert single-cell RNA sequencing data to Seurat-compatible format.

This script reads an H5AD file and converts it to the MTX format compatible with Seurat,
along with associated metadata files. All output files are gzipped.
"""

import os
import gzip
import shutil
import logging
from pathlib import Path
from typing import Union, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix

def get_project_root() -> Path:
    """Get the project root directory from the current script location.
    
    Returns:
        Path to the project root directory
    """
    current_file = Path(__file__).resolve()  # Get the path of the current script
    return current_file.parent.parent  # Go up two levels from preprocessing/ to main_directory/

def setup_logging(log_level: str = "INFO") -> None:
    """Configure logging for the script.

    Args:
        log_level: Desired logging level (default: "INFO")
    """
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def create_directory(directory: Union[str, Path]) -> Path:
    """Create directory if it doesn't exist.

    Args:
        directory: Path to create

    Returns:
        Path object of created directory
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    return directory

def read_h5ad_data(file_path: Union[str, Path]) -> ad.AnnData:
    """Read H5AD file and verify its contents.

    Args:
        file_path: Path to H5AD file

    Returns:
        AnnData object containing the single-cell data

    Raises:
        FileNotFoundError: If the input file doesn't exist
        ValueError: If the data is empty or invalid
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")

    adata = sc.read_h5ad(file_path)
    
    if adata.X.shape[0] == 0 or adata.X.shape[1] == 0:
        raise ValueError("The input data appears to be empty")
    
    return adata

def gzip_file(input_path: Path) -> None:
    """Gzip a file and remove the original.

    Args:
        input_path: Path to the file to be gzipped
    """
    output_path = Path(str(input_path) + '.gz')
    with open(input_path, 'rb') as f_in:
        with gzip.open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    input_path.unlink()  # Remove original file
    logging.info(f"Gzipped {input_path} to {output_path}")

def write_mtx_files(adata: ad.AnnData, output_dir: Path) -> None:
    """Write the count matrix and associated files in MTX format and gzip them.

    Args:
        adata: AnnData object containing the data
        output_dir: Directory to write the output files
    """
    mtx_dir = create_directory(output_dir / "mtx_files")
    
    # Write count matrix
    matrix_path = mtx_dir / "matrix.mtx"
    io.mmwrite(str(matrix_path), adata.X.T)
    gzip_file(matrix_path)
    
    # Write barcodes
    barcodes_path = mtx_dir / "barcodes.tsv"
    with open(barcodes_path, "w") as f:
        f.write('\n'.join(adata.obs_names) + '\n')
    gzip_file(barcodes_path)
    
    # Write features
    features_path = mtx_dir / "features.tsv"
    with open(features_path, "w") as f:
        f.write('\n'.join(adata.var_names) + '\n')
    gzip_file(features_path)

def write_metadata(adata: ad.AnnData, output_dir: Path) -> None:
    """Write observation and variable metadata to CSV files.

    Args:
        adata: AnnData object containing the data
        output_dir: Directory to write the output files
    """
    # Write metadata files
    metadata_path = output_dir / "metadata.csv"
    vardata_path = output_dir / "vardata.csv"
    
    adata.obs.to_csv(metadata_path)
    adata.var.to_csv(vardata_path)

def main() -> None:
    """Main function to execute the conversion workflow."""
    setup_logging()
    
    # Get project root directory
    project_root = get_project_root()
    
    # Define paths relative to project root
    input_file = project_root / "data/hlma_download/scRNA_h5ad/All_CellType_scsn_RNA.h5ad"
    output_dir = project_root / "data/processed/converted_files"
    
    try:
        # Read input data
        logging.info(f"Reading data from {input_file}")
        adata = read_h5ad_data(input_file)
        
        # Create output directory
        create_directory(output_dir)
        
        # Write MTX files
        logging.info("Writing and compressing MTX files")
        write_mtx_files(adata, output_dir)
        
        # Write metadata files
        logging.info("Writing metadata files")
        write_metadata(adata, output_dir)
        
        logging.info("Conversion completed successfully")
        
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise

if __name__ == "__main__":
    main()
