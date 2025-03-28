import sys
import os
import time
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                            QLabel, QPushButton, QFileDialog, QListWidget, QTabWidget, 
                            QComboBox, QGroupBox, QFormLayout, QLineEdit, QMessageBox,
                            QSplitter, QTreeWidget, QTreeWidgetItem, QCheckBox,
                            QAction, QMenuBar, QMenu, QDialog, QProgressBar, QSplashScreen)
from PyQt5.QtCore import Qt, QSize, QTimer
from PyQt5.QtGui import QColor, QFont, QPixmap, QIcon, QPainter

# For 3D visualization
from pyqtgraph.opengl import GLViewWidget, MeshData, GLMeshItem
import pyqtgraph.opengl as gl

# Biopython imports
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

class SplashScreen(QSplashScreen):
    def __init__(self):
        # Create pixmap for splash screen - average size
        splash_pix = QPixmap(650, 300)
        splash_pix.fill(Qt.white)
        super().__init__(splash_pix)
        
        # Create a separate widget for the progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFixedSize(400, 20)
        
        # Set the progress bar position
        self.progress_bar.move(50, 250)
        
        # Add content directly on the pixmap
        painter = QPainter(splash_pix)
        # Title
        painter.setFont(QFont("Arial", 22, QFont.Bold))
        painter.drawText(150, 120, "DeepMol v1.0")
        
        # Creator name
        painter.setFont(QFont("Arial", 9))
        painter.drawText(220, 200, "Deepanshu Sharma")
        painter.end()
        
        # Update the splash screen with the painted pixmap
        self.setPixmap(splash_pix)
        self.setWindowFlag(Qt.WindowStaysOnTopHint)
    
    def showProgress(self, value):
        """Show progress at the bottom of the splash screen"""
        self.showMessage(f"Loading... {value}%", Qt.AlignBottom | Qt.AlignCenter, Qt.black)

class CustomSelect(Select):
    """Custom selector for saving modified PDB structures"""
    def __init__(self, residues_to_remove=None, atoms_to_remove=None):
        self.residues_to_remove = residues_to_remove or set()
        self.atoms_to_remove = atoms_to_remove or set()
        
    def accept_residue(self, residue):
        residue_id = (residue.get_parent().id, residue.id[1])
        return residue_id not in self.residues_to_remove
    
    def accept_atom(self, atom):
        return atom.get_full_id() not in self.atoms_to_remove

class AboutDialog(QDialog):
    """Dialog showing information about the application"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About DeepMol")
        self.setFixedSize(400, 300)
        
        layout = QVBoxLayout()
        
        # App name and version
        title_label = QLabel("DeepMol v1.0")
        title_label.setFont(QFont("Arial", 18, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        
        # Creator info
        creator_label = QLabel("Deepanshu Sharma")
        creator_label.setFont(QFont("Arial", 12))
        creator_label.setAlignment(Qt.AlignCenter)
        
        # Description
        desc_label = QLabel("Molecular visualization and editing tool for PDB files.")
        desc_label.setWordWrap(True)
        desc_label.setAlignment(Qt.AlignCenter)
        
        # Features
        features_label = QLabel("Features:")
        features_label.setFont(QFont("Arial", 10, QFont.Bold))
        
        features_list = QLabel(
            "• Interactive 3D visualization of protein structures\n"
            "• Advanced structure editing capabilities\n"
            "• Intuitive atom and residue manipulation\n"
            "• Efficient 3D navigation and visualization options\n"
            "• PDB file format support"
        )
        
        # Close button
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        
        # Add all widgets to layout
        layout.addWidget(title_label)
        layout.addWidget(creator_label)
        layout.addSpacing(20)
        layout.addWidget(desc_label)
        layout.addSpacing(20)
        layout.addWidget(features_label)
        layout.addWidget(features_list)
        layout.addStretch(1)
        layout.addWidget(close_button)
        
        self.setLayout(layout)

class ProteinStructureExplorer(QMainWindow):
    def __init__(self):
        super().__init__()
        
        # Initialize variables
        self.structure = None
        self.current_chain = None
        self.current_model_id = 0
        self.deleted_residues = set()
        self.deleted_atoms = set()
        self.modified_coords = {}  # Store modified atom coordinates
        self.atom_meshes = []  # Store references to atom mesh items
        self.bond_lines = []   # Store references to bond line items
        
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle('DeepMol v1.0')
        self.setGeometry(100, 100, 1200, 800)
        
        # Create menu bar
        self.create_menu_bar()
        
        # Main widget and layout
        main_widget = QWidget()
        main_layout = QHBoxLayout(main_widget)
        
        # Create a splitter for resizable panels
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left panel - controls and structure info
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # File operations
        file_group = QGroupBox("File Operations")
        file_layout = QVBoxLayout()
        
        self.load_button = QPushButton("Load PDB File")
        self.load_button.clicked.connect(self.load_pdb)
        self.save_button = QPushButton("Save Modified Structure")
        self.save_button.clicked.connect(self.save_structure)
        self.save_button.setEnabled(False)
        
        file_layout.addWidget(self.load_button)
        file_layout.addWidget(self.save_button)
        file_group.setLayout(file_layout)
        left_layout.addWidget(file_group)
        
        # Structure Navigation
        nav_group = QGroupBox("Structure Navigation")
        nav_layout = QVBoxLayout()
        
        self.model_combo = QComboBox()
        self.model_combo.currentIndexChanged.connect(self.update_chain_list)
        
        self.chain_list = QListWidget()
        self.chain_list.itemClicked.connect(self.chain_selected)
        
        nav_layout.addWidget(QLabel("Model:"))
        nav_layout.addWidget(self.model_combo)
        nav_layout.addWidget(QLabel("Chains:"))
        nav_layout.addWidget(self.chain_list)
        nav_group.setLayout(nav_layout)
        left_layout.addWidget(nav_group)
        
        # Structure modification
        mod_group = QGroupBox("Structure Modification")
        mod_layout = QVBoxLayout()
        
        self.delete_residue_button = QPushButton("Delete Selected Residue")
        self.delete_residue_button.clicked.connect(self.delete_selected_residue)
        self.delete_residue_button.setEnabled(False)
        
        self.delete_atom_button = QPushButton("Delete Selected Atom")
        self.delete_atom_button.clicked.connect(self.delete_selected_atom)
        self.delete_atom_button.setEnabled(False)
        
        self.move_atom_button = QPushButton("Move Selected Atom")
        self.move_atom_button.clicked.connect(self.open_move_atom_dialog)
        self.move_atom_button.setEnabled(False)
        
        self.restore_button = QPushButton("Restore Original Structure")
        self.restore_button.clicked.connect(self.restore_structure)
        self.restore_button.setEnabled(False)
        
        mod_layout.addWidget(self.delete_residue_button)
        mod_layout.addWidget(self.delete_atom_button)
        mod_layout.addWidget(self.move_atom_button)
        mod_layout.addWidget(self.restore_button)
        mod_group.setLayout(mod_layout)
        left_layout.addWidget(mod_group)
        
        # Visualization Options
        viz_group = QGroupBox("Visualization Options")
        viz_layout = QVBoxLayout()
        
        self.show_backbone_only = QCheckBox("Show Backbone Only")
        self.show_backbone_only.toggled.connect(self.update_visualization)
        
        self.color_by_combo = QComboBox()
        self.color_by_combo.addItems(["Chain", "Residue Type", "Secondary Structure", "B-factor"])
        self.color_by_combo.currentIndexChanged.connect(self.update_visualization)
        
        viz_layout.addWidget(self.show_backbone_only)
        viz_layout.addWidget(QLabel("Color by:"))
        viz_layout.addWidget(self.color_by_combo)
        viz_group.setLayout(viz_layout)
        left_layout.addWidget(viz_group)
        
        # Middle panel - structure hierarchy and selection
        middle_panel = QWidget()
        middle_layout = QVBoxLayout(middle_panel)
        
        self.structure_tree = QTreeWidget()
        self.structure_tree.setHeaderLabels(["Structure Elements"])
        self.structure_tree.itemClicked.connect(self.element_selected)
        middle_layout.addWidget(QLabel("Structure Hierarchy:"))
        middle_layout.addWidget(self.structure_tree)
        
        # Right panel - 3D visualization
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        self.view_3d = GLViewWidget()
        self.view_3d.setCameraPosition(distance=50)
        right_layout.addWidget(self.view_3d)
        
        # Add all panels to splitter
        splitter.addWidget(left_panel)
        splitter.addWidget(middle_panel)
        splitter.addWidget(right_panel)
        
        # Set initial sizes
        splitter.setSizes([300, 300, 600])
        
        self.setCentralWidget(main_widget)
    
    def create_menu_bar(self):
        """Create the application menu bar"""
        menu_bar = self.menuBar()
        
        # File menu
        file_menu = menu_bar.addMenu('File')
        
        open_action = QAction('Open PDB', self)
        open_action.setShortcut('Ctrl+O')
        open_action.triggered.connect(self.load_pdb)
        
        save_action = QAction('Save', self)
        save_action.setShortcut('Ctrl+S')
        save_action.triggered.connect(self.save_structure)
        
        exit_action = QAction('Exit', self)
        exit_action.setShortcut('Alt+F4')
        exit_action.triggered.connect(self.close)
        
        file_menu.addAction(open_action)
        file_menu.addAction(save_action)
        file_menu.addSeparator()
        file_menu.addAction(exit_action)
        
        # Edit menu
        edit_menu = menu_bar.addMenu('Edit')
        
        restore_action = QAction('Restore Original', self)
        restore_action.triggered.connect(self.restore_structure)
        
        edit_menu.addAction(restore_action)
        
        # View menu
        view_menu = menu_bar.addMenu('View')
        
        backbone_action = QAction('Show Backbone Only', self)
        backbone_action.setCheckable(True)
        backbone_action.toggled.connect(self.toggle_backbone_view)
        
        view_menu.addAction(backbone_action)
        
        # Help menu
        help_menu = menu_bar.addMenu('Help')
        
        about_action = QAction('About', self)
        about_action.triggered.connect(self.show_about_dialog)
        
        help_menu.addAction(about_action)
    
    def toggle_backbone_view(self, state):
        """Toggle backbone-only view"""
        self.show_backbone_only.setChecked(state)
        self.update_visualization()
    
    def show_about_dialog(self):
        """Show the about dialog"""
        dialog = AboutDialog(self)
        dialog.exec_()
    
    def load_pdb(self):
        """Load a PDB file"""
        file_path, _ = QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB Files (*.pdb *.ent);;All Files (*)")
        
        if file_path:
            try:
                # Create a progress dialog
                progress = QProgressBar(self)
                progress.setRange(0, 100)
                progress.setValue(0)
                progress.setWindowTitle("Loading PDB File")
                progress.setGeometry(300, 300, 300, 40)
                progress.show()
                
                QApplication.processEvents()
                
                # Update progress
                progress.setValue(10)
                QApplication.processEvents()
                
                # Parse the structure
                parser = PDBParser(QUIET=True)
                progress.setValue(30)
                QApplication.processEvents()
                
                self.structure = parser.get_structure("protein", file_path)
                progress.setValue(60)
                QApplication.processEvents()
                
                # Store the original file path
                self.original_file_path = file_path
                
                # Clear previous data
                self.deleted_residues = set()
                self.deleted_atoms = set()
                self.modified_coords = {}
                
                # Update UI
                self.update_model_combo()
                self.update_structure_tree()
                self.save_button.setEnabled(True)
                self.restore_button.setEnabled(True)
                
                progress.setValue(90)
                QApplication.processEvents()
                
                # Display the structure
                self.current_model_id = 0
                self.update_visualization()
                
                # Enable modification buttons if structure is loaded
                if self.structure is not None:
                    self.delete_residue_button.setEnabled(True)
                    self.delete_atom_button.setEnabled(True)
                    self.move_atom_button.setEnabled(True)
                
                progress.setValue(100)
                QApplication.processEvents()
                progress.close()
                
                # Set window title to include filename
                filename = os.path.basename(file_path)
                self.setWindowTitle(f"DeepMol v1.0 - {filename}")
                
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Could not load PDB file: {str(e)}")
    
    def update_model_combo(self):
        """Update the model dropdown with available models"""
        self.model_combo.clear()
        if self.structure:
            for model in self.structure:
                self.model_combo.addItem(f"Model {model.id}")
    
    def update_chain_list(self):
        """Update the chain list based on the selected model"""
        self.chain_list.clear()
        if self.structure:
            model_id = self.model_combo.currentIndex()
            if model_id >= 0:
                self.current_model_id = model_id
                model = list(self.structure)[model_id]
                for chain in model:
                    self.chain_list.addItem(f"Chain {chain.id}")
                
                # Update visualization
                self.update_visualization()
    
    def chain_selected(self, item):
        """Handle chain selection"""
        if self.structure and item:
            chain_id = item.text().split()[-1]
            model = list(self.structure)[self.current_model_id]
            self.current_chain = model[chain_id]
            self.update_visualization(highlight_chain=chain_id)
    
    def update_structure_tree(self):
        """Update the structure hierarchy tree"""
        self.structure_tree.clear()
        if not self.structure:
            return
            
        for model in self.structure:
            model_item = QTreeWidgetItem(self.structure_tree)
            model_item.setText(0, f"Model {model.id}")
            model_item.setData(0, Qt.UserRole, ("model", model.id))
            
            for chain in model:
                chain_item = QTreeWidgetItem(model_item)
                chain_item.setText(0, f"Chain {chain.id}")
                chain_item.setData(0, Qt.UserRole, ("chain", model.id, chain.id))
                
                for residue in chain:
                    # Skip deleted residues
                    residue_id = (chain.id, residue.id[1])
                    if residue_id in self.deleted_residues:
                        continue
                        
                    res_name = residue.get_resname()
                    res_id = residue.id[1]
                    residue_item = QTreeWidgetItem(chain_item)
                    residue_item.setText(0, f"{res_name} {res_id}")
                    residue_item.setData(0, Qt.UserRole, ("residue", model.id, chain.id, residue.id))
                    
                    for atom in residue:
                        # Skip deleted atoms
                        if atom.get_full_id() in self.deleted_atoms:
                            continue
                            
                        atom_item = QTreeWidgetItem(residue_item)
                        atom_item.setText(0, f"{atom.name} {atom.element}")
                        atom_item.setData(0, Qt.UserRole, ("atom", model.id, chain.id, residue.id, atom.name))
        
        self.structure_tree.expandToDepth(1)
    
    def element_selected(self, item, column):
        """Handle selection in the structure tree"""
        if not item:
            return
            
        element_data = item.data(0, Qt.UserRole)
        if not element_data:
            return
            
        element_type = element_data[0]
        
        if element_type == "residue":
            self.delete_residue_button.setEnabled(True)
            self.delete_atom_button.setEnabled(False)
            self.move_atom_button.setEnabled(False)
        elif element_type == "atom":
            self.delete_residue_button.setEnabled(False)
            self.delete_atom_button.setEnabled(True)
            self.move_atom_button.setEnabled(True)
        else:
            self.delete_residue_button.setEnabled(False)
            self.delete_atom_button.setEnabled(False)
            self.move_atom_button.setEnabled(False)
            
        # Highlight the selected element
        self.highlight_selected_element(element_data)
    
    def highlight_selected_element(self, element_data):
        """Highlight the selected element in the 3D view"""
        element_type = element_data[0]
        
        # Implementation depends on the visualization approach
        # For now, just update the visualization with highlighting
        self.update_visualization(highlight_element=element_data)
    
    def delete_selected_residue(self):
        """Delete the selected residue"""
        item = self.structure_tree.currentItem()
        if not item:
            return
            
        element_data = item.data(0, Qt.UserRole)
        if not element_data or element_data[0] != "residue":
            return
            
        model_id, chain_id, residue_id = element_data[1:4]
        
        # Add to deleted residues set
        self.deleted_residues.add((chain_id, residue_id[1]))
        
        # Update UI
        self.update_structure_tree()
        self.update_visualization()
        
        QMessageBox.information(self, "Success", f"Residue {residue_id[1]} deleted. Changes will be applied when saving.")
    
    def delete_selected_atom(self):
        """Delete the selected atom"""
        item = self.structure_tree.currentItem()
        if not item:
            return
            
        element_data = item.data(0, Qt.UserRole)
        if not element_data or element_data[0] != "atom":
            return
            
        model_id, chain_id, residue_id, atom_name = element_data[1:5]
        
        # Get the atom's full ID
        model = list(self.structure)[model_id]
        chain = model[chain_id]
        residue = chain[residue_id]
        atom = residue[atom_name]
        atom_full_id = atom.get_full_id()
        
        # Add to deleted atoms set
        self.deleted_atoms.add(atom_full_id)
        
        # Update UI
        self.update_structure_tree()
        self.update_visualization()
        
        QMessageBox.information(self, "Success", f"Atom {atom_name} deleted. Changes will be applied when saving.")
    
    def open_move_atom_dialog(self):
        """Open dialog to move the selected atom"""
        item = self.structure_tree.currentItem()
        if not item:
            return
            
        element_data = item.data(0, Qt.UserRole)
        if not element_data or element_data[0] != "atom":
            return
            
        model_id, chain_id, residue_id, atom_name = element_data[1:5]
        
        # Get the atom
        model = list(self.structure)[model_id]
        chain = model[chain_id]
        residue = chain[residue_id]
        atom = residue[atom_name]
        
        # Create dialog
        dialog = QWidget(self, Qt.Window)
        dialog.setWindowTitle(f"Move Atom {atom_name}")
        dialog.setGeometry(300, 300, 400, 200)
        
        layout = QVBoxLayout(dialog)
        
        # Current coordinates
        current_coords = atom.get_coord()
        
        # Form for new coordinates
        form_layout = QFormLayout()
        
        x_input = QLineEdit(str(current_coords[0]))
        y_input = QLineEdit(str(current_coords[1]))
        z_input = QLineEdit(str(current_coords[2]))
        
        form_layout.addRow("X coordinate:", x_input)
        form_layout.addRow("Y coordinate:", y_input)
        form_layout.addRow("Z coordinate:", z_input)
        
        layout.addLayout(form_layout)
        
        # Apply button
        apply_button = QPushButton("Apply")
        
        def apply_move():
            try:
                new_x = float(x_input.text())
                new_y = float(y_input.text())
                new_z = float(z_input.text())
                
                new_coords = np.array([new_x, new_y, new_z], dtype=np.float32)
                
                # Store the modification
                atom_full_id = atom.get_full_id()
                self.modified_coords[atom_full_id] = new_coords
                
                # Apply the change
                atom.set_coord(new_coords)
                
                # Update visualization
                self.update_visualization()
                
                dialog.close()
                QMessageBox.information(self, "Success", f"Atom {atom_name} moved. Changes will be applied when saving.")
                
            except ValueError:
                QMessageBox.warning(dialog, "Error", "Please enter valid coordinates (numbers only).")
        
        apply_button.clicked.connect(apply_move)
        layout.addWidget(apply_button)
        
        dialog.show()
    
    def restore_structure(self):
        """Restore the original structure"""
        if not hasattr(self, 'original_file_path') or not self.original_file_path:
            QMessageBox.warning(self, "Error", "Original file path not found.")
            return
            
        # Create a fresh copy from the original file
        try:
            parser = PDBParser(QUIET=True)
            self.structure = parser.get_structure("protein", self.original_file_path)
            
            # Reset modifications
            self.deleted_residues = set()
            self.deleted_atoms = set()
            self.modified_coords = {}
            
            # Update UI
            self.update_model_combo()
            self.update_structure_tree()
            self.update_visualization()
            
            QMessageBox.information(self, "Success", "Structure restored to original state.")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to restore structure: {str(e)}")
    
    def save_structure(self):
        """Save the modified structure"""
        if not self.structure:
            return
            
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Modified PDB", "", "PDB Files (*.pdb);;All Files (*)")
        
        if file_path:
            try:
                # Apply atom coordinate modifications
                for atom_id, new_coords in self.modified_coords.items():
                    # Navigate to the atom using its full ID
                    model_id, chain_id, residue_id, atom_name = atom_id[1], atom_id[2], atom_id[3], atom_id[4][0]
                    model = self.structure[model_id]
                    chain = model[chain_id]
                    residue = chain[residue_id]
                    atom = residue[atom_name]
                    
                    # Set the new coordinates
                    atom.set_coord(new_coords)
                
                # Create a selector that excludes deleted residues and atoms
                selector = CustomSelect(
                    residues_to_remove=self.deleted_residues,
                    atoms_to_remove=self.deleted_atoms
                )
                
                # Save using the selector
                io = PDBIO()
                io.set_structure(self.structure)
                io.save(file_path, selector)
                
                QMessageBox.information(self, "Success", f"Modified structure saved to {file_path}")
                
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Failed to save structure: {str(e)}")
    
    def update_visualization(self, highlight_chain=None, highlight_element=None):
        """Update the 3D visualization of the structure"""
        # Clear current visualization
        self.view_3d.clear()
        
        # Clear stored references
        self.atom_meshes = []
        self.bond_lines = []
        
        if not self.structure:
            return
            
        # Get current visualization options
        backbone_only = self.show_backbone_only.isChecked()
        color_by = self.color_by_combo.currentText()
        
        # Get the current model
        model = list(self.structure)[self.current_model_id]
        
        # Define colors for chains
        chain_colors = {
            'A': (1.0, 0.0, 0.0, 1.0),  # Red
            'B': (0.0, 0.0, 1.0, 1.0),  # Blue
            'C': (0.0, 1.0, 0.0, 1.0),  # Green
            'D': (1.0, 1.0, 0.0, 1.0),  # Yellow
            'E': (1.0, 0.0, 1.0, 1.0),  # Magenta
            'F': (0.0, 1.0, 1.0, 1.0),  # Cyan
            'G': (0.5, 0.5, 0.5, 1.0),  # Gray
            'H': (1.0, 0.5, 0.0, 1.0)   # Orange
        }
        
        # Define colors for residue types
        residue_colors = {
            # Hydrophobic
            'ALA': (0.8, 0.8, 0.8, 1.0),  # Light gray
            'VAL': (0.8, 0.8, 0.8, 1.0),
            'LEU': (0.8, 0.8, 0.8, 1.0),
            'ILE': (0.8, 0.8, 0.8, 1.0),
            'MET': (0.8, 0.8, 0.8, 1.0),
            'PHE': (0.8, 0.8, 0.8, 1.0),
            'TRP': (0.8, 0.8, 0.8, 1.0),
            # Polar
            'SER': (0.0, 1.0, 1.0, 1.0),  # Cyan
            'THR': (0.0, 1.0, 1.0, 1.0),
            'ASN': (0.0, 1.0, 1.0, 1.0),
            'GLN': (0.0, 1.0, 1.0, 1.0),
            'TYR': (0.0, 1.0, 1.0, 1.0),
            # Basic
            'LYS': (0.0, 0.0, 1.0, 1.0),  # Blue
            'ARG': (0.0, 0.0, 1.0, 1.0),
            'HIS': (0.0, 0.0, 1.0, 1.0),
            # Acidic
            'ASP': (1.0, 0.0, 0.0, 1.0),  # Red
            'GLU': (1.0, 0.0, 0.0, 1.0),
            # Special
            'CYS': (1.0, 1.0, 0.0, 1.0),  # Yellow
            'PRO': (0.5, 0.5, 0.0, 1.0),  # Olive
            'GLY': (1.0, 0.5, 0.5, 1.0)   # Pink
        }
        
        # Optimization: Limit the number of atoms displayed
        max_atoms_per_chain = 1000  # Adjust this value based on performance
        
        # Process each chain and collect data for visualization
        atom_data = []
        bond_data = []
        backbone_atoms = ['N', 'CA', 'C', 'O']
        
        for chain in model:
            # Skip if we're highlighting a specific chain and this isn't it
            if highlight_chain and chain.id != highlight_chain:
                continue
                
            chain_atom_count = 0
            prev_atom = None
            
            for residue in chain:
                # Skip deleted residues
                residue_id = (chain.id, residue.id[1])
                if residue_id in self.deleted_residues:
                    continue
                
                # Determine residue color based on selected coloring scheme
                if color_by == "Chain":
                    residue_color = chain_colors.get(chain.id, (0.5, 0.5, 0.5, 1.0))
                elif color_by == "Residue Type":
                    residue_color = residue_colors.get(residue.get_resname(), (0.5, 0.5, 0.5, 1.0))
                else:
                    residue_color = (0.5, 0.5, 0.5, 1.0)
                
                for atom in residue:
                    # Skip deleted atoms
                    if atom.get_full_id() in self.deleted_atoms:
                        continue
                    
                    # Skip non-backbone atoms if backbone_only is checked
                    if backbone_only and atom.name not in backbone_atoms:
                        continue
                    
                    # Get atom position
                    pos = atom.get_coord()
                    
                    # Check if this atom has modified coordinates
                    atom_full_id = atom.get_full_id()
                    if atom_full_id in self.modified_coords:
                        pos = self.modified_coords[atom_full_id]
                    
                    # Ensure pos is a numpy array
                    pos = np.array(pos, dtype=np.float32)
                    
                    # Check if this atom is highlighted
                    is_highlighted = False
                    if highlight_element:
                        element_type = highlight_element[0]
                        if element_type == "atom":
                            # Check if this is the highlighted atom
                            _, model_id, chain_id, res_id, atom_name = highlight_element
                            if (model.id == model_id and 
                                chain.id == chain_id and 
                                residue.id == res_id and 
                                atom.name == atom_name):
                                is_highlighted = True
                        elif element_type == "residue":
                            # Check if this atom belongs to the highlighted residue
                            _, model_id, chain_id, res_id = highlight_element
                            if (model.id == model_id and 
                                chain.id == chain_id and 
                                residue.id == res_id):
                                is_highlighted = True
                    
                    # Set atom color
                    if is_highlighted:
                        atom_color = (1.0, 1.0, 1.0, 1.0)  # White for highlighted atoms
                    else:
                        # Use element-specific colors for atoms
                        element = atom.element
                        if element == 'C':
                            atom_color = (0.2, 0.2, 0.2, 1.0)
                        elif element == 'N':
                            atom_color = (0.0, 0.0, 1.0, 1.0)
                        elif element == 'O':
                            atom_color = (1.0, 0.0, 0.0, 1.0)
                        elif element == 'S':
                            atom_color = (1.0, 1.0, 0.0, 1.0)
                        else:
                            atom_color = residue_color
                    
                    # Determine atom size based on element
                    if atom.element == 'H':
                        size = 0.5
                    elif atom.element in ['C', 'N', 'O']:
                        size = 1.0
                    else:
                        size = 1.2
                    
                    # Make highlighted atoms larger
                    if is_highlighted:
                        size *= 1.5
                    
                    # Store atom data for rendering
                    atom_data.append({
                        'pos': pos,
                        'color': atom_color,
                        'size': size
                    })
                    
                    # Check for backbone bonds
                    if prev_atom is not None and atom.name in backbone_atoms and prev_atom.name in backbone_atoms:
                        # Only connect sequential backbone atoms
                        if ((prev_atom.name == 'N' and atom.name == 'CA') or
                            (prev_atom.name == 'CA' and atom.name == 'C') or
                            (prev_atom.name == 'C' and atom.name == 'O')):
                            
                            prev_pos = prev_atom.get_coord()
                            
                            # Check if previous atom has modified coordinates
                            prev_full_id = prev_atom.get_full_id()
                            if prev_full_id in self.modified_coords:
                                prev_pos = self.modified_coords[prev_full_id]
                            
                            # Ensure prev_pos is a numpy array
                            prev_pos = np.array(prev_pos, dtype=np.float32)
                            
                            # Add bond
                            bond_data.append({
                                'pos': np.array([prev_pos, pos]),
                                'color': residue_color
                            })
                    
                    # Update previous atom
                    if atom.name in backbone_atoms:
                        prev_atom = atom
                    
                    chain_atom_count += 1
                    
                # Limit atoms per chain for performance
                if chain_atom_count >= max_atoms_per_chain:
                    break
        
        # Create atom spheres
        level_of_detail = 8  # Lower for faster rendering
        for data in atom_data:
            sphere = gl.GLMeshItem(
                meshdata=gl.MeshData.sphere(rows=level_of_detail, cols=level_of_detail, radius=data['size']),
                smooth=True,
                color=data['color'],
                shader='shaded',
                glOptions='opaque'
            )
            sphere.translate(data['pos'][0], data['pos'][1], data['pos'][2])
            self.view_3d.addItem(sphere)
            self.atom_meshes.append(sphere)
        
        # Create bonds
        for bond in bond_data:
            line = gl.GLLinePlotItem(
                pos=bond['pos'],
                color=bond['color'],
                width=2,
                antialias=True
            )
            self.view_3d.addItem(line)
            self.bond_lines.append(line)
        
        # Add coordinate axes for reference
        axis_length = 10
        axes = [
            ((0, 0, 0), (axis_length, 0, 0), (1, 0, 0, 1)),  # X-axis (red)
            ((0, 0, 0), (0, axis_length, 0), (0, 1, 0, 1)),  # Y-axis (green)
            ((0, 0, 0), (0, 0, axis_length), (0, 0, 1, 1))   # Z-axis (blue)
        ]
        
        for start, end, color in axes:
            line = gl.GLLinePlotItem(
                pos=np.array([start, end]),
                color=color,
                width=2,
                antialias=True
            )
            self.view_3d.addItem(line)
            self.bond_lines.append(line)

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("DeepMol")
    
    # Create splash screen
    splash = SplashScreen()
    splash.show()
    app.processEvents()
    
    # Create main window but don't show yet
    window = ProteinStructureExplorer()
    
    # Simple delay with progress updates
    for i in range(0, 101, 10):
        splash.showProgress(i)
        app.processEvents()
        time.sleep(0.1)  # Short delay
    
    # Short delay before showing main window
    time.sleep(0.5)
    
    # Show main window
    window.show()
    splash.finish(window)
    
    # Execute application
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()