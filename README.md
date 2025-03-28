# DeepMol v1.0

Molecular visualization and editing tool for PDB files

## Features

- **Interactive 3D Visualization**: View protein structures in a fully interactive 3D environment
- **Structure Editing**: Add, delete, or modify atoms and residues with precision
- **Intuitive Navigation**: Easily explore complex protein structures through a hierarchical view
- **Multiple Visualization Options**: View proteins by chain, residue type, and more
- **PDB Support**: Load and save standard PDB files with all modifications preserved

## Installation from Source

```bash
# Clone the repository
git clone https://github.com/your-username/deepmol.git
cd deepmol

# Install dependencies
pip install -r requirements.txt

# Run the application
python protein_structure_explorer.py
```

## Dependencies

- Python 3.7+
- PyQt5
- PyQtGraph
- NumPy
- Biopython

## Usage

1. **Loading a Structure**: Use File → Open or the "Load PDB File" button to open a PDB file
2. **Navigation**: 
   - Use the tree view to navigate the hierarchical structure
   - Select chains from the dropdown menu
   - Rotate the 3D view by clicking and dragging
   - Zoom with the scroll wheel
3. **Editing**:
   - Select an atom or residue in the tree view
   - Use the modification buttons to delete or move elements
   - Save your changes with File → Save
4. **Visualization Options**:
   - Toggle "Show Backbone Only" for a clearer view
   - Select different coloring schemes from the dropdown

## Building from Source

To create a standalone executable:

```bash
# Install PyInstaller
pip install pyinstaller

# Create executable
pyinstaller --onefile --windowed DeepMol.py
```

The executable will be created in the `dist` directory.

## Author

Deepanshu Sharma

## Acknowledgments

- PyMOL and Chimera for inspiration
- [Biopython](https://biopython.org/) for PDB parsing capabilities
- The entire open-source community for their invaluable tools and libraries
