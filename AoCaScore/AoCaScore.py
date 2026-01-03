"""
Aortic Valve Calcium Score - Guided Workflow Plugin

A 3D Slicer extension for quantification of aortic valve calcification using
the Agatston scoring method. Features a guided tab-based workflow with automated
segmentation, classification, and comprehensive PDF reporting.

Author: Vittorio Censullo
Institution: AITeRTC (Associazione Italiana Tecnici di Radiologia Esperti in TC)
Year: 2025
Version: 3.0 (Non-Contrast CT Only)
License: MIT

Based on ACC/AHA guidelines and JAHA 2024 literature:

STANDARD AGATSTON SCORING (Non-Contrast CT):
- Fixed threshold: 130 HU (standard for non-contrast cardiac CT)
- Sex-specific severity thresholds (JAHA 2024):
  * Women: Severe ≥1300 AU, Moderate 400-1300 AU
  * Men: Severe ≥2000 AU, Moderate 1000-2000 AU
- 2D slice-by-slice labeling (matches Syngo.via and commercial software)
- Automatic handling of overlapping slices (3mm thickness/1.5mm increment)
- Minimum lesion criterion: 2 contiguous pixels per slice (Agatston standard)
- Density-based weighting (factors 1-4)
- 3D visualization and PDF reporting

Key Features:
- Validated against Syngo.via (±1.3% accuracy)
- Automatic slice overlap detection and correction
- Per-slice density factor calculation
- Comprehensive PDF reporting with patient information
- ROI-based segmentation for precise calcium detection

Repository: https://github.com/vcensullo/SlicerCaScore
"""

import os
import unittest
import vtk
import qt
import ctk
import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
import numpy as np
from datetime import datetime

#
# SlicerCaScore
#

class AoCaScore(ScriptedLoadableModule):
    """Aortic Valve Calcium Score Module"""

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Aortic Valve Calcium Score"
        self.parent.categories = ["Cardiac", "Quantification"]
        self.parent.dependencies = []
        self.parent.contributors = ["Vittorio Censullo (AITeRTC - 2025)"]
        self.parent.helpText = """
        <h3>Aortic Valve Calcium Score - Guided Workflow (v2.0)</h3>

        <p>Quantification of aortic valve calcification using the Agatston scoring method with support for both non-contrast and contrast-enhanced CT.</p>

        <h4>Features:</h4>
        <ul>
        <li>Guided tab-based workflow for easy use</li>
        <li><b>Dual-mode support:</b> Non-contrast CT (fixed threshold) and Contrast-Enhanced CT (dynamic threshold)</li>
        <li>Automatic threshold calculation for contrast CT (Mean Aorta HU + 4×SD)</li>
        <li>ROI-based or manual paint segmentation</li>
        <li>Automated Agatston score calculation with volume filtering</li>
        <li>Sex-specific severity classification adapted for each CT modality</li>
        <li>3D visualization with density-based coloring</li>
        <li>Comprehensive PDF reporting with screenshots</li>
        </ul>

        <h4>CT Mode Selection:</h4>
        <ul>
        <li><b>Non-Contrast CT:</b> Standard Agatston (130 HU threshold, women ≥1300 AU / men ≥2000 AU for severe AS)</li>
        <li><b>Contrast-Enhanced CT:</b> Modified Agatston (dynamic threshold, women ≥1430 AU / men ≥1840 AU for severe AS)</li>
        </ul>

        <h4>Author Information:</h4>
        <ul>
        <li><b>Author:</b> Vittorio Censullo</li>
        <li><b>Institution:</b> AITeRTC (Artificial Intelligence in Teleradiology and Telemedicine Research Center)</li>
        <li><b>Year:</b> 2025</li>
        <li><b>Version:</b> 2.0 (Dual-Mode)</li>
        <li><b>License:</b> MIT</li>
        <li><b>Repository:</b> <a href="https://github.com/vcensullo/AoCaScore">GitHub</a></li>
        </ul>

        <h4>References:</h4>
        <p>Based on ACC/AHA guidelines and recent literature:</p>
        <ul>
        <li>Non-contrast thresholds: JAHA 2024</li>
        <li>Contrast thresholds: Circulation: Cardiovascular Imaging 2025</li>
        <li>Dynamic threshold method: J Clin Med 2024, EuroIntervention</li>
        </ul>
        """
        self.parent.acknowledgementText = """
        <p>Developed by <b>Vittorio Censullo</b> at <b>AITeRTC</b> (Associazione Italiana Tecnici di Radiologia Esperti in TC) - 2025</p>

        <p>This plugin provides standardized aortic valve calcium scoring according to
        ACC/AHA guidelines and recent literature, including sex-specific severity thresholds.</p>

        <p>For issues, suggestions, or contributions, please visit the
        <a href="https://github.com/vcensullo/AoCaScore">GitHub repository</a>.</p>
        """

#
# SlicerCaScoreWidget
#

class AoCaScoreWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Main widget with tab-based guided workflow"""

    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)

        self.logic = None
        self.volumeNode = None
        self.segmentationNode = None
        self.roiNode = None
        self.roiTransformNode = None

        # Default threshold (can be modified in Settings tab)
        self.thresholdValue = 130  # HU for non-contrast CT (Standard Agatston)

        # Patient data
        self.patientInfo = {
            'name': '',
            'id': '',
            'sex': 'M',
            'date': datetime.now().strftime("%Y-%m-%d"),
            'age': ''
        }

    def setup(self):
        """Setup UI"""
        ScriptedLoadableModuleWidget.setup(self)

        self.logic = AoCaScoreLogic()

        # Module Header - Simple
        headerFrame = qt.QFrame()
        headerFrame.setStyleSheet("""
            QFrame {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 #0078d7, stop:1 #005a9e);
                border-radius: 8px;
                padding: 15px;
                margin-bottom: 10px;
            }
        """)
        headerLayout = qt.QVBoxLayout(headerFrame)
        self.layout.addWidget(headerFrame)

        titleLabel = qt.QLabel("💙 Aortic Valve Calcium Score")
        titleLabel.setStyleSheet("color: white; font-size: 20px; font-weight: bold;")
        headerLayout.addWidget(titleLabel)

        # Patient info section (collapsible, always visible)
        self.setupPatientInfoSection()

        # Create tab widget for workflow steps
        self.tabWidget = qt.QTabWidget()
        self.layout.addWidget(self.tabWidget)

        # Create tabs
        self.setupSettingsTab()  # Settings & Dependencies
        self.setupStep1Tab()  # Setup & View
        self.setupStep2Tab()  # Segmentation
        self.setupStep3Tab()  # Results
        self.setupStep4Tab()  # Report

        # Disable tabs initially (enable as we progress)
        # Settings tab is always enabled
        self.tabWidget.setTabEnabled(2, False)  # Segmentation
        self.tabWidget.setTabEnabled(3, False)  # Results
        self.tabWidget.setTabEnabled(4, False)  # Report

        self.layout.addStretch(1)

        # Add observer to auto-select volumes when they are added to the scene
        self.sceneObserverTag = slicer.mrmlScene.AddObserver(
            slicer.mrmlScene.NodeAddedEvent,
            self.onNodeAdded
        )

        # Check if there are already volumes in the scene and auto-select the first one
        # Use QTimer to delay this slightly to ensure UI is fully initialized
        qt.QTimer.singleShot(100, self.checkAndSelectVolume)

    def setupPatientInfoSection(self):
        """Patient information section"""
        patientBox = ctk.ctkCollapsibleButton()
        patientBox.text = "Patient Information"
        patientBox.collapsed = True
        self.layout.addWidget(patientBox)
        patientLayout = qt.QFormLayout(patientBox)

        # Patient Name
        self.patientNameEdit = qt.QLineEdit()
        self.patientNameEdit.setPlaceholderText("Enter patient name")
        patientLayout.addRow("Name:", self.patientNameEdit)

        # Patient ID
        self.patientIDEdit = qt.QLineEdit()
        self.patientIDEdit.setPlaceholderText("Enter patient ID")
        patientLayout.addRow("ID:", self.patientIDEdit)

        # Sex
        sexWidget = qt.QWidget()
        sexLayout = qt.QHBoxLayout(sexWidget)
        sexLayout.setContentsMargins(0, 0, 0, 0)
        self.sexButtonGroup = qt.QButtonGroup()
        self.maleButton = qt.QRadioButton("Male")
        self.femaleButton = qt.QRadioButton("Female")
        self.maleButton.setChecked(True)
        self.sexButtonGroup.addButton(self.maleButton, 0)
        self.sexButtonGroup.addButton(self.femaleButton, 1)
        sexLayout.addWidget(self.maleButton)
        sexLayout.addWidget(self.femaleButton)
        patientLayout.addRow("Sex:", sexWidget)

        # Age
        self.patientAgeEdit = qt.QLineEdit()
        self.patientAgeEdit.setPlaceholderText("Age")
        patientLayout.addRow("Age:", self.patientAgeEdit)

    def setupSettingsTab(self):
        """Tab 0: Settings - Dependencies and Configuration"""
        settingsTab = qt.QWidget()
        layout = qt.QVBoxLayout(settingsTab)

        # === Section: Dependencies ===
        depsFrame = qt.QFrame()
        depsFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        depsLayout = qt.QVBoxLayout(depsFrame)
        layout.addWidget(depsFrame)

        depsTitle = qt.QLabel("<b>Dependencies Management</b>")
        depsTitle.setStyleSheet("font-size: 14px; color: #0078d7;")
        depsLayout.addWidget(depsTitle)

        depsInfo = qt.QLabel(
            "This plugin requires the 'reportlab' library for PDF report generation.\n"
            "Click the button below to install missing dependencies."
        )
        depsInfo.setWordWrap(True)
        depsLayout.addWidget(depsInfo)

        # Check dependencies status
        self.depsStatusLabel = qt.QLabel()
        self.depsStatusLabel.setWordWrap(True)
        depsLayout.addWidget(self.depsStatusLabel)

        # Install button (create BEFORE calling updateDependenciesStatus)
        installButtonLayout = qt.QHBoxLayout()
        self.installDepsButton = qt.QPushButton("Install/Update Dependencies")
        self.installDepsButton.setStyleSheet("QPushButton { background-color: #0078d7; color: white; padding: 8px; font-weight: bold; }")
        self.installDepsButton.clicked.connect(self.onInstallDependencies)
        installButtonLayout.addWidget(self.installDepsButton)
        installButtonLayout.addStretch()
        depsLayout.addLayout(installButtonLayout)

        # Now check dependencies status
        self.updateDependenciesStatus()

        depsLayout.addWidget(qt.QLabel(""))  # Spacer

        # === Section: Threshold Configuration ===
        thresholdFrame = qt.QFrame()
        thresholdFrame.setStyleSheet("QFrame { background-color: #fff3cd; border-radius: 5px; padding: 10px; margin: 5px; }")
        thresholdLayout = qt.QVBoxLayout(thresholdFrame)
        layout.addWidget(thresholdFrame)

        thresholdTitle = qt.QLabel("<b>Threshold Configuration</b>")
        thresholdTitle.setStyleSheet("font-size: 14px; color: #856404;")
        thresholdLayout.addWidget(thresholdTitle)

        thresholdInfo = qt.QLabel(
            "<b>WARNING:</b> The default threshold of 130 HU is the standard for non-contrast CT calcium detection.<br>"
            "Modifying this value should only be done with proper clinical justification."
        )
        thresholdInfo.setWordWrap(True)
        thresholdInfo.setStyleSheet("color: #856404;")
        thresholdLayout.addWidget(thresholdInfo)

        # Threshold slider
        thresholdControlLayout = qt.QFormLayout()
        self.thresholdSlider = ctk.ctkSliderWidget()
        self.thresholdSlider.singleStep = 1
        self.thresholdSlider.minimum = 100
        self.thresholdSlider.maximum = 200
        self.thresholdSlider.value = self.thresholdValue
        self.thresholdSlider.suffix = " HU"
        self.thresholdSlider.decimals = 0
        self.thresholdSlider.connect('valueChanged(double)', self.onThresholdChanged)
        thresholdControlLayout.addRow("Minimum Threshold:", self.thresholdSlider)
        thresholdLayout.addLayout(thresholdControlLayout)

        # Current value display
        self.currentThresholdLabel = qt.QLabel(f"Current threshold: <b>{self.thresholdValue} HU</b>")
        self.currentThresholdLabel.setStyleSheet("color: #856404; font-size: 12px;")
        thresholdLayout.addWidget(self.currentThresholdLabel)

        # Reset button
        resetButtonLayout = qt.QHBoxLayout()
        resetButton = qt.QPushButton("Reset to Default (130 HU)")
        resetButton.clicked.connect(self.onResetThreshold)
        resetButtonLayout.addWidget(resetButton)
        resetButtonLayout.addStretch()
        thresholdLayout.addLayout(resetButtonLayout)

        # === Section: Plugin Information ===
        layout.addWidget(qt.QLabel(""))  # Spacer

        infoFrame = qt.QFrame()
        infoFrame.setStyleSheet("QFrame { background-color: #e7f3ff; border-radius: 5px; padding: 10px; margin: 5px; }")
        infoLayout = qt.QVBoxLayout(infoFrame)
        layout.addWidget(infoFrame)

        infoTitle = qt.QLabel("<b>Plugin Information</b>")
        infoTitle.setStyleSheet("font-size: 14px; color: #0078d7;")
        infoLayout.addWidget(infoTitle)

        infoText = qt.QLabel(
            "<b>Aortic Valve Calcium Score Plugin</b><br><br>"
            "Version: 3.0 (Non-Contrast Only - Syngo.via Validated)<br>"
            "Author: Vittorio Censullo<br>"
            "Institution: AITeRTC<br>"
            "Year: 2025<br><br>"
            "Purpose: Quantification of aortic valve calcification using Standard Agatston scoring<br><br>"
            "<b>CT Protocol:</b><br>"
            "• Non-Contrast ECG-gated CT<br>"
            "• 3mm slice thickness (1.5mm increment supported)<br>"
            "• 130 HU threshold (adjustable 100-200 HU)<br>"
            "• 2D slice-by-slice labeling (Agatston standard)<br><br>"
            "<b>Workflow:</b><br>"
            "1. Setup View: Load CT volume and position viewing planes<br>"
            "2. Calcium Segmentation: Place ROI and segment calcium automatically<br>"
            "3. Results & Analysis: Calculate Agatston score (validated vs Syngo.via)<br>"
            "4. Generate Report: Create PDF report with all findings<br><br>"
            "<b>Severity Thresholds (JAHA 2024):</b><br>"
            "Women: Severe ≥1300 AU, Moderate 400-1300 AU<br>"
            "Men: Severe ≥2000 AU, Moderate 1000-2000 AU"
        )
        infoText.setWordWrap(True)
        infoLayout.addWidget(infoText)

        layout.addStretch()

        # Add tab
        self.tabWidget.addTab(settingsTab, "⚙ Settings")

    def setupStep1Tab(self):
        """Tab 1: Setup View - Slice intersections for positioning"""
        tab1 = qt.QWidget()
        layout = qt.QVBoxLayout(tab1)

        # === Section: Select CT Volume ===
        volumeFrame = qt.QFrame()
        volumeFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        volumeLayout = qt.QVBoxLayout(volumeFrame)
        layout.addWidget(volumeFrame)

        volumeTitle = qt.QLabel("<b>Select CT Volume</b>")
        volumeTitle.setStyleSheet("font-size: 14px; color: #0078d7;")
        volumeLayout.addWidget(volumeTitle)

        # CT Volume selector
        selectorFrame = qt.QFrame()
        selectorLayout = qt.QFormLayout(selectorFrame)
        volumeLayout.addWidget(selectorFrame)

        self.volumeSelector = slicer.qMRMLNodeComboBox()
        self.volumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.volumeSelector.selectNodeUponCreation = False
        self.volumeSelector.addEnabled = False
        self.volumeSelector.removeEnabled = False
        self.volumeSelector.noneEnabled = True
        self.volumeSelector.showHidden = False
        self.volumeSelector.setMRMLScene(slicer.mrmlScene)
        self.volumeSelector.setToolTip("Select the cardiac CT volume")
        selectorLayout.addRow("CT Volume:", self.volumeSelector)

        # Volume info
        self.volumeInfoLabel = qt.QLabel("No volume selected")
        self.volumeInfoLabel.setWordWrap(True)
        self.volumeInfoLabel.setStyleSheet("QLabel { padding: 5px; background-color: #ffffff; border: 1px solid #dee2e6; border-radius: 3px; }")
        selectorLayout.addRow("Info:", self.volumeInfoLabel)

        # === Section: Enable Slice Intersections ===
        sliceFrame = qt.QFrame()
        sliceFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        sliceLayout = qt.QVBoxLayout(sliceFrame)
        layout.addWidget(sliceFrame)

        sliceTitle = qt.QLabel("<b>Slice Intersection Planes</b>")
        sliceTitle.setStyleSheet("font-size: 14px; color: #0078d7;")
        sliceLayout.addWidget(sliceTitle)

        instructionLabel = qt.QLabel(
            "<b>Instructions:</b><br>"
            "1. Enable slice intersections to visualize plane orientations<br>"
            "2. Click on the aortic valve in any slice view<br>"
            "3. The intersection planes will help you position accurately<br>"
            "4. Apply the optimal layout when ready"
        )
        instructionLabel.setWordWrap(True)
        instructionLabel.setStyleSheet("QLabel { padding: 10px; background-color: #e3f2fd; border-left: 4px solid #2196f3; border-radius: 3px; }")
        sliceLayout.addWidget(instructionLabel)

        # === Section: Apply Layout ===
        layoutFrame = qt.QFrame()
        layoutFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        layoutLayout = qt.QVBoxLayout(layoutFrame)
        layout.addWidget(layoutFrame)

        layoutTitle = qt.QLabel("<b>Apply Optimal Layout</b>")
        layoutTitle.setStyleSheet("font-size: 14px; color: #0078d7;")
        layoutLayout.addWidget(layoutTitle)

        # Window/Level option
        self.autoWindowLevelCheckbox = qt.QCheckBox("Apply calcium preset window/level (W:1500, L:300)")
        self.autoWindowLevelCheckbox.setChecked(True)
        self.autoWindowLevelCheckbox.setStyleSheet("font-weight: bold;")
        layoutLayout.addWidget(self.autoWindowLevelCheckbox)

        # Apply layout button
        self.setupViewButton = qt.QPushButton("Apply Optimal Layout")
        self.setupViewButton.setIcon(qt.QIcon(":/Icons/LayoutFourUpView.png"))
        self.setupViewButton.setIconSize(qt.QSize(32, 32))
        self.setupViewButton.setMinimumHeight(50)
        self.setupViewButton.enabled = False
        self.setupViewButton.setStyleSheet("""
            QPushButton {
                background-color: #28a745;
                color: white;
                font-weight: bold;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #218838;
            }
            QPushButton:disabled {
                background-color: #6c757d;
            }
        """)
        layoutLayout.addWidget(self.setupViewButton)

        # === Status Section ===
        statusFrame = qt.QFrame()
        statusFrame.setStyleSheet("QFrame { background-color: #e8f5e9; border: 2px solid #4caf50; border-radius: 5px; padding: 10px; margin: 5px; }")
        statusLayout = qt.QVBoxLayout(statusFrame)
        layout.addWidget(statusFrame)

        self.step1StatusLabel = qt.QLabel("Status: Waiting for volume selection")
        self.step1StatusLabel.setStyleSheet("color: #2e7d32; font-weight: bold;")
        statusLayout.addWidget(self.step1StatusLabel)

        layout.addStretch(1)

        # Add tab
        self.tabWidget.addTab(tab1, "1. Setup View")

        # Connect signals
        self.volumeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onVolumeSelected)
        self.setupViewButton.connect('clicked(bool)', self.onApplyLayout)

    def setupStep2Tab(self):
        """Tab 2: Calcium Segmentation"""
        tab2 = qt.QWidget()
        layout = qt.QVBoxLayout(tab2)

        # === Section Title ===
        sectionTitle = qt.QLabel("<b>Choose Segmentation Method:</b>")
        sectionTitle.setStyleSheet("font-size: 14px; color: #333; margin-top: 15px; margin-bottom: 5px;")
        layout.addWidget(sectionTitle)

        # === Method Selection Buttons (3 large icon buttons) ===
        methodsFrame = qt.QFrame()
        methodsLayout = qt.QHBoxLayout(methodsFrame)
        methodsLayout.setSpacing(10)
        layout.addWidget(methodsFrame)

        # Method 1: ROI Button
        self.roiMethodButton = qt.QPushButton()
        self.roiMethodButton.setCheckable(True)
        self.roiMethodButton.setMinimumSize(150, 120)
        self.roiMethodButton.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                border: 2px solid #cccccc;
                border-radius: 10px;
                padding: 15px;
            }
            QPushButton:hover {
                border: 2px solid #0078d7;
                background-color: #f0f8ff;
            }
            QPushButton:checked {
                border: 3px solid #0078d7;
                background-color: #e3f2fd;
            }
            QPushButton:disabled {
                background-color: #f5f5f5;
                border: 2px solid #e0e0e0;
            }
        """)
        self.roiMethodButton.enabled = False
        roiLayout = qt.QVBoxLayout()
        roiIcon = qt.QLabel()
        # Load custom icon from module's Icons folder
        modulePath = os.path.dirname(slicer.util.modulePath(self.__module__))
        iconPath = os.path.join(modulePath, "Resources", "Icons", "roi_box.png")
        if os.path.exists(iconPath):
            roiIcon.setPixmap(qt.QPixmap(iconPath).scaled(48, 48, qt.Qt.KeepAspectRatio, qt.Qt.SmoothTransformation))
        else:
            roiIcon.setText("📦")  # Fallback emoji
            roiIcon.setStyleSheet("font-size: 36px;")
        roiIcon.setAlignment(qt.Qt.AlignCenter)
        roiText = qt.QLabel("<b>ROI Box</b><br><small>Automatic threshold<br>in region</small>")
        roiText.setAlignment(qt.Qt.AlignCenter)
        roiText.setStyleSheet("color: #333;")
        roiLayout.addWidget(roiIcon)
        roiLayout.addWidget(roiText)
        self.roiMethodButton.setLayout(roiLayout)
        methodsLayout.addWidget(self.roiMethodButton)

        # Method 2: Click & Grow Button
        self.clickGrowMethodButton = qt.QPushButton()
        self.clickGrowMethodButton.setCheckable(True)
        self.clickGrowMethodButton.setMinimumSize(150, 120)
        self.clickGrowMethodButton.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                border: 2px solid #28a745;
                border-radius: 10px;
                padding: 15px;
            }
            QPushButton:hover {
                border: 2px solid #28a745;
                background-color: #f0fff4;
            }
            QPushButton:checked {
                border: 3px solid #28a745;
                background-color: #e8f5e9;
            }
            QPushButton:disabled {
                background-color: #f5f5f5;
                border: 2px solid #e0e0e0;
            }
        """)
        self.clickGrowMethodButton.enabled = False
        clickGrowLayout = qt.QVBoxLayout()
        clickGrowIcon = qt.QLabel()
        # Load custom icon from module's Icons folder
        iconPath = os.path.join(modulePath, "Resources", "Icons", "click_grow.png")
        if os.path.exists(iconPath):
            clickGrowIcon.setPixmap(qt.QPixmap(iconPath).scaled(48, 48, qt.Qt.KeepAspectRatio, qt.Qt.SmoothTransformation))
        else:
            clickGrowIcon.setText("🖱️")  # Fallback emoji
            clickGrowIcon.setStyleSheet("font-size: 36px;")
        clickGrowIcon.setAlignment(qt.Qt.AlignCenter)
        clickGrowText = qt.QLabel("<b>Click & Grow</b><br><small>⭐ Recommended<br>One-click segmentation</small>")
        clickGrowText.setAlignment(qt.Qt.AlignCenter)
        clickGrowText.setStyleSheet("color: #28a745; font-weight: bold;")
        clickGrowLayout.addWidget(clickGrowIcon)
        clickGrowLayout.addWidget(clickGrowText)
        self.clickGrowMethodButton.setLayout(clickGrowLayout)
        methodsLayout.addWidget(self.clickGrowMethodButton)

        # Method 3: Brush Button
        self.brushMethodButton = qt.QPushButton()
        self.brushMethodButton.setCheckable(True)
        self.brushMethodButton.setMinimumSize(150, 120)
        self.brushMethodButton.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                border: 2px solid #cccccc;
                border-radius: 10px;
                padding: 15px;
            }
            QPushButton:hover {
                border: 2px solid #0078d7;
                background-color: #f0f8ff;
            }
            QPushButton:checked {
                border: 3px solid #0078d7;
                background-color: #e3f2fd;
            }
            QPushButton:disabled {
                background-color: #f5f5f5;
                border: 2px solid #e0e0e0;
            }
        """)
        self.brushMethodButton.enabled = False
        brushLayout = qt.QVBoxLayout()
        brushIcon = qt.QLabel()
        # Load custom icon from module's Icons folder
        iconPath = os.path.join(modulePath, "Resources", "Icons", "brush.png")
        if os.path.exists(iconPath):
            brushIcon.setPixmap(qt.QPixmap(iconPath).scaled(48, 48, qt.Qt.KeepAspectRatio, qt.Qt.SmoothTransformation))
        else:
            brushIcon.setText("🖌️")  # Fallback emoji
            brushIcon.setStyleSheet("font-size: 36px;")
        brushIcon.setAlignment(qt.Qt.AlignCenter)
        brushText = qt.QLabel("<b>Brush</b><br><small>Manual paint<br>& erase</small>")
        brushText.setAlignment(qt.Qt.AlignCenter)
        brushText.setStyleSheet("color: #333;")
        brushLayout.addWidget(brushIcon)
        brushLayout.addWidget(brushText)
        self.brushMethodButton.setLayout(brushLayout)
        methodsLayout.addWidget(self.brushMethodButton)

        # === Method 1: ROI Box + Threshold (Collapsible) ===
        self.method1Collapsible = ctk.ctkCollapsibleButton()
        self.method1Collapsible.text = "ROI Box Method - Controls"
        self.method1Collapsible.collapsed = True
        self.method1Collapsible.visible = False
        layout.addWidget(self.method1Collapsible)

        method1Layout = qt.QVBoxLayout(self.method1Collapsible)

        method1Label = qt.QLabel("<b>Method 1: ROI Box + Automatic Threshold</b>")
        method1Label.setStyleSheet("font-size: 13px; color: #0078d7;")
        method1Layout.addWidget(method1Label)

        method1Desc = qt.QLabel("Define a box around the aortic valve, then automatically segment all voxels ≥ 130 HU (fixed threshold for non-contrast CT)")
        method1Desc.setWordWrap(True)
        method1Desc.setStyleSheet("color: #666666;")
        method1Layout.addWidget(method1Desc)

        roiButtonsFrame = qt.QFrame()
        roiButtonsLayout = qt.QHBoxLayout(roiButtonsFrame)
        method1Layout.addWidget(roiButtonsFrame)

        self.placeROIButton = qt.QPushButton("Place ROI Box")
        self.placeROIButton.setCheckable(True)
        self.placeROIButton.setIcon(qt.QIcon(":/Icons/MarkupsROIPlace.png"))
        self.placeROIButton.setIconSize(qt.QSize(20, 20))
        self.placeROIButton.enabled = False
        roiButtonsLayout.addWidget(self.placeROIButton)

        self.applyROIThresholdButton = qt.QPushButton("Apply Threshold in ROI")
        self.applyROIThresholdButton.setIcon(qt.QIcon(":/Icons/Go.png"))
        self.applyROIThresholdButton.setIconSize(qt.QSize(20, 20))
        self.applyROIThresholdButton.enabled = False
        roiButtonsLayout.addWidget(self.applyROIThresholdButton)

        # Show ROI checkbox
        self.showROICheckbox = qt.QCheckBox("Show ROI")
        self.showROICheckbox.setChecked(True)
        self.showROICheckbox.enabled = False
        method1Layout.addWidget(self.showROICheckbox)

        # ROI Rotation Controls
        rotationLabel = qt.QLabel("<i>ROI Rotation (degrees):</i>")
        rotationLabel.setStyleSheet("color: #666666; margin-top: 5px;")
        method1Layout.addWidget(rotationLabel)

        rotationFrame = qt.QFrame()
        rotationLayout = qt.QFormLayout(rotationFrame)
        method1Layout.addWidget(rotationFrame)

        self.roiRotationLR = ctk.ctkSliderWidget()
        self.roiRotationLR.singleStep = 1
        self.roiRotationLR.minimum = -180
        self.roiRotationLR.maximum = 180
        self.roiRotationLR.value = 0
        self.roiRotationLR.enabled = False
        self.roiRotationLR.setToolTip("Rotate ROI around Left-Right axis")
        rotationLayout.addRow("L-R:", self.roiRotationLR)

        self.roiRotationPA = ctk.ctkSliderWidget()
        self.roiRotationPA.singleStep = 1
        self.roiRotationPA.minimum = -180
        self.roiRotationPA.maximum = 180
        self.roiRotationPA.value = 0
        self.roiRotationPA.enabled = False
        self.roiRotationPA.setToolTip("Rotate ROI around Posterior-Anterior axis")
        rotationLayout.addRow("P-A:", self.roiRotationPA)

        self.roiRotationIS = ctk.ctkSliderWidget()
        self.roiRotationIS.singleStep = 1
        self.roiRotationIS.minimum = -180
        self.roiRotationIS.maximum = 180
        self.roiRotationIS.value = 0
        self.roiRotationIS.enabled = False
        self.roiRotationIS.setToolTip("Rotate ROI around Inferior-Superior axis")
        rotationLayout.addRow("I-S:", self.roiRotationIS)

        self.resetRotationButton = qt.QPushButton("Reset Rotation")
        self.resetRotationButton.setIcon(qt.QIcon(":/Icons/Undo.png"))
        self.resetRotationButton.enabled = False
        method1Layout.addWidget(self.resetRotationButton)

        # === Method 2: Click and Grow (One-Click Segmentation) - RECOMMENDED ===
        self.method2Collapsible = ctk.ctkCollapsibleButton()
        self.method2Collapsible.text = "Click & Grow Method - Controls"
        self.method2Collapsible.collapsed = True
        self.method2Collapsible.visible = False
        layout.addWidget(self.method2Collapsible)

        method2Layout = qt.QVBoxLayout(self.method2Collapsible)

        method2Label = qt.QLabel("<b>✨ Method 2: Click & Grow (RECOMMENDED)</b>")
        method2Label.setStyleSheet("font-size: 14px; color: #28a745; font-weight: bold;")
        method2Layout.addWidget(method2Label)

        method2Desc = qt.QLabel("⚡ Quick and easy! Click on each calcification to segment it automatically with 3D propagation.<br><i>Similar to Syngo.via workflow</i>")
        method2Desc.setWordWrap(True)
        method2Desc.setStyleSheet("color: #666666; margin: 5px 0;")
        method2Layout.addWidget(method2Desc)

        self.clickGrowButton = qt.QPushButton("  Enable Click & Grow")
        self.clickGrowButton.setCheckable(True)
        self.clickGrowButton.setIcon(qt.QIcon(":/Icons/MarkupsFiducialMouseModePlace.png"))
        self.clickGrowButton.setIconSize(qt.QSize(28, 28))
        self.clickGrowButton.setMinimumHeight(45)
        self.clickGrowButton.enabled = False
        self.clickGrowButton.setToolTip("Click on calcification to segment it automatically")
        self.clickGrowButton.setStyleSheet("""
            QPushButton {
                background-color: #28a745;
                color: white;
                font-size: 13px;
                font-weight: bold;
                border-radius: 5px;
                padding: 10px;
                text-align: left;
            }
            QPushButton:hover {
                background-color: #218838;
            }
            QPushButton:checked {
                background-color: #c82333;
            }
            QPushButton:disabled {
                background-color: #cccccc;
                color: #666666;
            }
        """)
        method2Layout.addWidget(self.clickGrowButton)

        # === Method 3: Manual Paint ===
        self.method3Collapsible = ctk.ctkCollapsibleButton()
        self.method3Collapsible.text = "Brush Method - Controls"
        self.method3Collapsible.collapsed = True
        self.method3Collapsible.visible = False
        layout.addWidget(self.method3Collapsible)

        method3Layout = qt.QVBoxLayout(self.method3Collapsible)

        method3Label = qt.QLabel("<b>Method 3: Manual Paint/Erase</b>")
        method3Label.setStyleSheet("font-size: 13px; color: #0078d7;")
        method3Layout.addWidget(method3Label)

        method3Desc = qt.QLabel("Manually paint or erase calcifications with brush, only voxels ≥ 130 HU will be segmented")
        method3Desc.setWordWrap(True)
        method3Desc.setStyleSheet("color: #666666;")
        method3Layout.addWidget(method3Desc)

        paintButtonsFrame = qt.QFrame()
        paintButtonsLayout = qt.QHBoxLayout(paintButtonsFrame)
        method3Layout.addWidget(paintButtonsFrame)

        self.paintButton = qt.QPushButton("Paint Mode")
        self.paintButton.setCheckable(True)
        self.paintButton.setIcon(qt.QIcon(":/Icons/SegmentEditorPaint.png"))
        self.paintButton.setIconSize(qt.QSize(20, 20))
        self.paintButton.enabled = False
        paintButtonsLayout.addWidget(self.paintButton)

        self.eraseButton = qt.QPushButton("Erase Mode")
        self.eraseButton.setCheckable(True)
        self.eraseButton.setIcon(qt.QIcon(":/Icons/SegmentEditorErase.png"))
        self.eraseButton.setIconSize(qt.QSize(20, 20))
        self.eraseButton.enabled = False
        paintButtonsLayout.addWidget(self.eraseButton)

        # === Parameters ===
        paramsFrame = qt.QFrame()
        paramsFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        paramsLayout = qt.QFormLayout(paramsFrame)
        layout.addWidget(paramsFrame)

        # Threshold info (displays current threshold from Settings)
        self.thresholdInfoLabel = qt.QLabel(f"{self.thresholdValue} HU")
        self.thresholdInfoLabel.setStyleSheet("font-weight: bold; color: #0078d7;")
        paramsLayout.addRow("Threshold:", self.thresholdInfoLabel)

        thresholdNote = qt.QLabel("(Modifiable in Settings tab)")
        thresholdNote.setStyleSheet("font-size: 9px; color: #666;")
        paramsLayout.addRow("", thresholdNote)

        # Brush size slider
        self.brushSizeSlider = ctk.ctkSliderWidget()
        self.brushSizeSlider.singleStep = 1
        self.brushSizeSlider.minimum = 1
        self.brushSizeSlider.maximum = 20
        self.brushSizeSlider.value = 5
        self.brushSizeSlider.suffix = " mm"
        self.brushSizeSlider.setToolTip("Brush diameter for painting")
        paramsLayout.addRow("Brush Size:", self.brushSizeSlider)

        # === Editing Tools ===
        editFrame = qt.QFrame()
        editLayout = qt.QHBoxLayout(editFrame)
        layout.addWidget(editFrame)

        self.clearSegmentButton = qt.QPushButton("Clear All")
        self.clearSegmentButton.enabled = False
        editLayout.addWidget(self.clearSegmentButton)

        self.undoButton = qt.QPushButton("Undo")
        self.undoButton.enabled = False
        editLayout.addWidget(self.undoButton)

        # Show 3D preview
        self.preview3DCheckbox = qt.QCheckBox("Show 3D Preview")
        self.preview3DCheckbox.setChecked(True)
        layout.addWidget(self.preview3DCheckbox)

        layout.addStretch(1)

        # Add tab
        self.tabWidget.addTab(tab2, "2. Calcium Segmentation")

        # Connect signals
        # Method selection buttons
        self.roiMethodButton.connect('clicked(bool)', self.onROIMethodSelected)
        self.clickGrowMethodButton.connect('clicked(bool)', self.onClickGrowMethodSelected)
        self.brushMethodButton.connect('clicked(bool)', self.onBrushMethodSelected)

        self.placeROIButton.connect('clicked(bool)', self.onPlaceROIToggled)
        self.applyROIThresholdButton.connect('clicked(bool)', self.onApplyROIThreshold)
        self.showROICheckbox.connect('toggled(bool)', self.onToggleROIVisibility)
        self.roiRotationLR.connect('valueChanged(double)', self.onROIRotationChanged)
        self.roiRotationPA.connect('valueChanged(double)', self.onROIRotationChanged)
        self.roiRotationIS.connect('valueChanged(double)', self.onROIRotationChanged)
        self.resetRotationButton.connect('clicked(bool)', self.onResetROIRotation)
        self.clickGrowButton.connect('clicked(bool)', self.onClickGrowToggled)
        self.paintButton.connect('clicked(bool)', self.onPaintToggled)
        self.eraseButton.connect('clicked(bool)', self.onEraseToggled)
        self.clearSegmentButton.connect('clicked(bool)', self.onClearSegmentation)
        self.preview3DCheckbox.connect('toggled(bool)', self.onPreview3DToggled)

    # ============================================================
    # Settings Tab Callbacks
    # ============================================================

    def updateDependenciesStatus(self):
        """Check and display dependencies status"""
        # Check both reportlab and matplotlib
        reportlabInstalled = False
        matplotlibInstalled = False

        try:
            import reportlab
            reportlabVersion = reportlab.Version
            reportlabInstalled = True
        except ImportError:
            reportlabVersion = None

        try:
            import matplotlib
            matplotlibVersion = matplotlib.__version__
            matplotlibInstalled = True
        except ImportError:
            matplotlibVersion = None

        # Build status message
        status = []
        if reportlabInstalled:
            status.append(f"✓ reportlab {reportlabVersion}")
        else:
            status.append("✗ reportlab NOT installed")

        if matplotlibInstalled:
            status.append(f"✓ matplotlib {matplotlibVersion}")
        else:
            status.append("✗ matplotlib NOT installed")

        statusText = "\n".join(status)

        # Update UI
        if reportlabInstalled and matplotlibInstalled:
            self.depsStatusLabel.setText(statusText)
            self.depsStatusLabel.setStyleSheet("color: green; font-weight: bold;")
            self.installDepsButton.setText("Update Dependencies")
        elif reportlabInstalled or matplotlibInstalled:
            statusText += "\n\n⚠ Some dependencies are missing"
            self.depsStatusLabel.setText(statusText)
            self.depsStatusLabel.setStyleSheet("color: orange; font-weight: bold;")
            self.installDepsButton.setText("Install Missing Dependencies")
        else:
            statusText += "\n\n✗ PDF reports and Bull's Eye charts will not work"
            self.depsStatusLabel.setText(statusText)
            self.depsStatusLabel.setStyleSheet("color: red; font-weight: bold;")
            self.installDepsButton.setText("Install Dependencies")

    def onInstallDependencies(self):
        """Install/update required dependencies"""
        try:
            import subprocess
            import sys

            # Show progress message
            progressDialog = qt.QProgressDialog("Installing dependencies...", "Cancel", 0, 0)
            progressDialog.setWindowModality(qt.Qt.WindowModal)
            progressDialog.setCancelButton(None)
            progressDialog.setLabelText("Installing reportlab and matplotlib...")
            progressDialog.show()
            slicer.app.processEvents()

            # Get Python executable from Slicer
            pythonExe = sys.executable

            # Install both reportlab and matplotlib
            packages = ["reportlab", "matplotlib"]
            results = []
            errors = []

            for package in packages:
                progressDialog.setLabelText(f"Installing {package}...")
                slicer.app.processEvents()

                result = subprocess.run(
                    [pythonExe, "-m", "pip", "install", "--upgrade", package],
                    capture_output=True,
                    text=True
                )

                if result.returncode == 0:
                    results.append(f"✓ {package} installed successfully")
                else:
                    errors.append(f"✗ {package}: {result.stderr}")

            progressDialog.close()

            # Update status
            self.updateDependenciesStatus()

            # Show results
            if not errors:
                slicer.util.infoDisplay(
                    "All dependencies installed successfully!\n\n" +
                    "\n".join(results) +
                    "\n\nYou can now:\n" +
                    "• Generate PDF reports (reportlab)\n" +
                    "• Create Bull's Eye charts (matplotlib)"
                )
            else:
                slicer.util.warningDisplay(
                    "Installation completed with some errors:\n\n" +
                    "\n".join(results) + "\n\n" +
                    "\n".join(errors) + "\n\n" +
                    "Please try installing manually using:\n" +
                    f"{pythonExe} -m pip install reportlab matplotlib"
                )

        except Exception as e:
            if 'progressDialog' in locals():
                progressDialog.close()
            slicer.util.errorDisplay(f"Error during installation:\n{str(e)}")

    def onThresholdChanged(self, value):
        """Update threshold value"""
        self.thresholdValue = int(value)
        self.currentThresholdLabel.setText(f"Current threshold: <b>{self.thresholdValue} HU</b>")
        # Update threshold display in Segmentation tab if it exists
        if hasattr(self, 'thresholdInfoLabel'):
            self.thresholdInfoLabel.setText(f"{self.thresholdValue} HU")

    def onResetThreshold(self):
        """Reset threshold to default value"""
        self.thresholdValue = 130
        self.thresholdSlider.value = 130
        self.updateThresholdDisplay()
        # Update threshold display in Segmentation tab if it exists
        if hasattr(self, 'thresholdInfoLabel'):
            self.thresholdInfoLabel.setText(f"{self.thresholdValue} HU")
        slicer.util.infoDisplay("Threshold reset to default: 130 HU")

    # ============================================================
    # Step 3 Tab Setup
    # ============================================================

    def setupStep3Tab(self):
        """Tab 3: Results & Analysis"""
        tab3 = qt.QWidget()
        layout = qt.QVBoxLayout(tab3)

        # === Calculate Button ===
        calcFrame = qt.QFrame()
        calcFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        calcLayout = qt.QVBoxLayout(calcFrame)
        layout.addWidget(calcFrame)

        self.calculateButton = qt.QPushButton("Calculate Agatston Score")
        self.calculateButton.setIcon(qt.QIcon(":/Icons/MarkupsGeneric.png"))
        self.calculateButton.setIconSize(qt.QSize(32, 32))
        self.calculateButton.setMinimumHeight(50)
        self.calculateButton.enabled = False
        self.calculateButton.setStyleSheet("""
            QPushButton {
                background-color: #28a745;
                color: white;
                font-weight: bold;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #218838;
            }
            QPushButton:disabled {
                background-color: #6c757d;
            }
        """)
        calcLayout.addWidget(self.calculateButton)

        # === Results Display ===
        resultsGroup = qt.QGroupBox("Quantitative Results")
        resultsLayout = qt.QVBoxLayout(resultsGroup)
        layout.addWidget(resultsGroup)

        # Main score display
        self.agatstonScoreLabel = qt.QLabel("Agatston Score: ---")
        self.agatstonScoreLabel.setStyleSheet("""
            QLabel {
                font-size: 24px;
                font-weight: bold;
                padding: 15px;
                background-color: #e8f4f8;
                border: 2px solid #0078d7;
                border-radius: 5px;
            }
        """)
        self.agatstonScoreLabel.setAlignment(qt.Qt.AlignCenter)
        resultsLayout.addWidget(self.agatstonScoreLabel)

        # Classification
        self.classificationLabel = qt.QLabel("Classification: ---")
        self.classificationLabel.setStyleSheet("""
            QLabel {
                font-size: 16px;
                font-weight: bold;
                padding: 10px;
                border-radius: 5px;
            }
        """)
        self.classificationLabel.setAlignment(qt.Qt.AlignCenter)
        resultsLayout.addWidget(self.classificationLabel)

        # Detailed metrics table
        self.resultsTable = qt.QTableWidget()
        self.resultsTable.setColumnCount(2)
        self.resultsTable.setHorizontalHeaderLabels(["Metric", "Value"])
        self.resultsTable.horizontalHeader().setStretchLastSection(True)
        self.resultsTable.setMinimumHeight(200)
        resultsLayout.addWidget(self.resultsTable)

        # === Visualization Options ===
        vizGroup = qt.QGroupBox("Visualization Options")
        vizLayout = qt.QVBoxLayout(vizGroup)
        layout.addWidget(vizGroup)

        vizButtonsFrame = qt.QFrame()
        vizButtonsLayout = qt.QHBoxLayout(vizButtonsFrame)
        vizLayout.addWidget(vizButtonsFrame)

        self.show3DButton = qt.QPushButton("Show 3D Model")
        self.show3DButton.setIcon(qt.QIcon(":/Icons/View3D.png"))
        self.show3DButton.setIconSize(qt.QSize(20, 20))
        self.show3DButton.enabled = False
        vizButtonsLayout.addWidget(self.show3DButton)

        self.showChartsButton = qt.QPushButton("Show Charts")
        self.showChartsButton.setIcon(qt.QIcon(":/Icons/chart.png"))
        self.showChartsButton.setIconSize(qt.QSize(20, 20))
        self.showChartsButton.enabled = False
        vizButtonsLayout.addWidget(self.showChartsButton)

        self.showBullsEyeButton = qt.QPushButton("Show Bull's Eye")
        self.showBullsEyeButton.setIcon(qt.QIcon(":/Icons/chart.png"))
        self.showBullsEyeButton.setIconSize(qt.QSize(20, 20))
        self.showBullsEyeButton.enabled = False
        vizButtonsLayout.addWidget(self.showBullsEyeButton)

        # 3D model options
        self.colorByDensityCheckbox = qt.QCheckBox("Color by density (Agatston weighting)")
        self.colorByDensityCheckbox.setChecked(True)
        vizLayout.addWidget(self.colorByDensityCheckbox)

        layout.addStretch(1)

        # Add tab
        self.tabWidget.addTab(tab3, "3. Results & Analysis")

        # Connect signals
        self.calculateButton.connect('clicked(bool)', self.onCalculate)
        self.show3DButton.connect('clicked(bool)', self.onShow3D)
        self.showChartsButton.connect('clicked(bool)', self.onShowCharts)
        self.showBullsEyeButton.connect('clicked(bool)', self.onShowBullsEye)

    def setupStep4Tab(self):
        """Tab 4: Generate Report"""
        tab4 = qt.QWidget()
        layout = qt.QVBoxLayout(tab4)

        # === Report Preview ===
        previewFrame = qt.QFrame()
        previewFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        previewLayout = qt.QVBoxLayout(previewFrame)
        layout.addWidget(previewFrame)

        previewLabel = qt.QLabel("<b>Report will include:</b>")
        previewLabel.setStyleSheet("font-size: 13px; color: #0078d7;")
        previewLayout.addWidget(previewLabel)

        reportItemsList = qt.QListWidget()
        reportItems = [
            "Patient information",
            "Acquisition parameters",
            "Screenshots (Axial, Sagittal, 3D view)",
            "Quantitative results table",
            "Agatston score and classification",
            "Density distribution chart",
            "3D visualization of calcifications",
            "Literature-based interpretation"
        ]
        for item in reportItems:
            reportItemsList.addItem(item)
        reportItemsList.setMaximumHeight(200)
        previewLayout.addWidget(reportItemsList)

        # === Report Options ===
        optionsGroup = qt.QGroupBox("Report Options")
        optionsLayout = qt.QFormLayout(optionsGroup)
        layout.addWidget(optionsGroup)

        self.includeScreenshotsCheckbox = qt.QCheckBox()
        self.includeScreenshotsCheckbox.setChecked(True)
        optionsLayout.addRow("Include Screenshots:", self.includeScreenshotsCheckbox)

        self.include3DCheckbox = qt.QCheckBox()
        self.include3DCheckbox.setChecked(True)
        optionsLayout.addRow("Include 3D Model:", self.include3DCheckbox)

        self.includeChartsCheckbox = qt.QCheckBox()
        self.includeChartsCheckbox.setChecked(True)
        optionsLayout.addRow("Include Charts:", self.includeChartsCheckbox)

        # === Output Path ===
        outputFrame = qt.QFrame()
        outputLayout = qt.QHBoxLayout(outputFrame)
        layout.addWidget(outputFrame)

        outputLabel = qt.QLabel("Output Path:")
        outputLayout.addWidget(outputLabel)

        self.outputPathEdit = qt.QLineEdit()
        self.outputPathEdit.setPlaceholderText("Select output directory...")
        outputLayout.addWidget(self.outputPathEdit)

        self.browseButton = qt.QPushButton("Browse")
        self.browseButton.setIcon(qt.QIcon(":/Icons/AnnotationROI.png"))
        self.browseButton.setIconSize(qt.QSize(20, 20))
        outputLayout.addWidget(self.browseButton)

        # === Generate Button ===
        generateFrame = qt.QFrame()
        generateFrame.setStyleSheet("QFrame { background-color: #f8f9fa; border-radius: 5px; padding: 10px; margin: 5px; }")
        generateLayout = qt.QVBoxLayout(generateFrame)
        layout.addWidget(generateFrame)

        # Info label
        self.reportStatusLabel = qt.QLabel("Status: Calculate Agatston Score first (Tab 3)")
        self.reportStatusLabel.setStyleSheet("color: #6c757d; font-weight: bold; padding: 5px;")
        generateLayout.addWidget(self.reportStatusLabel)

        self.generateReportButton = qt.QPushButton("Generate PDF Report")
        self.generateReportButton.setIcon(qt.QIcon(":/Icons/AnnotationROI.png"))
        self.generateReportButton.setIconSize(qt.QSize(32, 32))
        self.generateReportButton.setMinimumHeight(50)
        self.generateReportButton.enabled = False
        self.generateReportButton.setStyleSheet("""
            QPushButton {
                background-color: #28a745;
                color: white;
                font-weight: bold;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #218838;
            }
            QPushButton:disabled {
                background-color: #6c757d;
            }
        """)
        generateLayout.addWidget(self.generateReportButton)

        layout.addStretch(1)

        # Add tab
        self.tabWidget.addTab(tab4, "4. Generate Report")

        # Connect signals
        self.browseButton.connect('clicked(bool)', self.onBrowseOutput)
        self.generateReportButton.connect('clicked(bool)', self.onGenerateReport)

    # ========== CALLBACK METHODS ==========

    def onVolumeSelected(self):
        """Handle volume selection"""
        self.volumeNode = self.volumeSelector.currentNode()

        if self.volumeNode:
            # Update info
            spacing = self.volumeNode.GetSpacing()
            dimensions = self.volumeNode.GetImageData().GetDimensions()
            info = f"Dimensions: {dimensions[0]}x{dimensions[1]}x{dimensions[2]}\n"
            info += f"Spacing: {spacing[0]:.2f}x{spacing[1]:.2f}x{spacing[2]:.2f} mm"
            self.volumeInfoLabel.text = info

            # Enable setup button
            self.setupViewButton.enabled = True

            # Enable Segmentation tab immediately when volume is selected
            self.tabWidget.setTabEnabled(2, True)

            # Create segmentation automatically or reuse existing one
            if not self.segmentationNode:
                # Check if there's already a segmentation node in the scene
                try:
                    existingNode = slicer.util.getNode("AorticValveCalcium")
                    self.segmentationNode = existingNode
                    print("✓ Reusing existing segmentation node")
                except slicer.util.MRMLNodeNotFoundException:
                    self.segmentationNode = self.logic.createSegmentation(self.volumeNode)
                    print("✓ Segmentation initialized automatically")

            # Enable method selection buttons
            self.roiMethodButton.enabled = True
            self.clickGrowMethodButton.enabled = True
            self.brushMethodButton.enabled = True

            # Enable segmentation tools
            self.placeROIButton.enabled = True
            self.clickGrowButton.enabled = True
            self.paintButton.enabled = True
            self.eraseButton.enabled = True
            self.clearSegmentButton.enabled = True

            # Update status
            self.step1StatusLabel.setText("Status: Volume selected - Ready for segmentation (Tab 2)")
            self.step1StatusLabel.setStyleSheet("color: #28a745; font-weight: bold;")
        else:
            self.volumeInfoLabel.text = "No volume selected"
            self.setupViewButton.enabled = False

            # Disable Segmentation tab if no volume
            self.tabWidget.setTabEnabled(2, False)

            # Disable method buttons
            self.roiMethodButton.enabled = False
            self.clickGrowMethodButton.enabled = False
            self.brushMethodButton.enabled = False

            self.step1StatusLabel.setText("Status: Waiting for volume selection")
            self.step1StatusLabel.setStyleSheet("color: #6c757d; font-weight: bold;")

    def onNodeAdded(self, caller, event, calldata=None):
        """Called when a new node is added to the scene"""
        # Check if it's a volume node
        if calldata and calldata.IsA('vtkMRMLScalarVolumeNode'):
            # Auto-select this volume if no volume is currently selected
            if not self.volumeSelector.currentNode():
                self.volumeSelector.setCurrentNode(calldata)
                print(f"✓ Auto-selected newly added volume: {calldata.GetName()}")
                # Manually trigger the volume selection callback
                self.onVolumeSelected()

    def checkAndSelectVolume(self):
        """Check if there are volumes in the scene and auto-select if needed"""
        # Only auto-select if no volume is currently selected
        currentNode = self.volumeSelector.currentNode()
        if not currentNode:
            volumeNodes = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')
            if volumeNodes:
                # Select the most recently added volume (last in list)
                # Force the selector to update by blocking signals temporarily
                self.volumeSelector.blockSignals(True)
                self.volumeSelector.setCurrentNode(volumeNodes[-1])
                self.volumeSelector.blockSignals(False)
                print(f"✓ Auto-selected volume: {volumeNodes[-1].GetName()}")
                # Manually trigger the volume selection callback
                self.onVolumeSelected()
            else:
                print("ℹ No volumes found in scene")

    def enter(self):
        """Called when the user enters this module"""
        # Check and auto-select volume when user switches to this module
        qt.QTimer.singleShot(50, self.checkAndSelectVolume)

    def cleanup(self):
        """Cleanup when module is closed"""
        # Remove scene observer
        if hasattr(self, 'sceneObserverTag') and self.sceneObserverTag:
            slicer.mrmlScene.RemoveObserver(self.sceneObserverTag)
            self.sceneObserverTag = None

    def onApplyLayout(self):
        """Apply optimal layout and enable segmentation tab"""
        if not self.volumeNode:
            slicer.util.errorDisplay("Please select a CT volume first")
            return

        # Update patient info
        self.updatePatientInfo()

        # Apply layout
        self.logic.setupOptimalView(
            self.volumeNode,
            applyWindowLevel=self.autoWindowLevelCheckbox.isChecked()
        )

        # Update status
        self.step1StatusLabel.setText("Status: View configured - Ready for segmentation!")
        self.step1StatusLabel.setStyleSheet("color: #28a745; font-weight: bold;")

        slicer.util.infoDisplay(
            "View configured successfully!\n\n"
            "Layout: Four-Up View (Axial, Sagittal, Coronal, 3D)\n"
            + ("Window/Level applied: W=1500, L=300\n" if self.autoWindowLevelCheckbox.isChecked() else "") +
            "\nYou can now proceed to calcium segmentation (Tab 2)."
        )

    # ============================================================
    # Method Selection Callbacks
    # ============================================================

    def onROIMethodSelected(self, checked):
        """Toggle ROI method panel"""
        if checked:
            # Show ROI method panel, hide others
            self.method1Collapsible.visible = True
            self.method1Collapsible.collapsed = False
            self.method2Collapsible.visible = False
            self.method3Collapsible.visible = False
            # Uncheck other buttons
            self.clickGrowMethodButton.setChecked(False)
            self.brushMethodButton.setChecked(False)
        else:
            # Hide ROI method panel
            self.method1Collapsible.visible = False

    def onClickGrowMethodSelected(self, checked):
        """Toggle Click & Grow method panel"""
        if checked:
            # Show Click & Grow method panel, hide others
            self.method1Collapsible.visible = False
            self.method2Collapsible.visible = True
            self.method2Collapsible.collapsed = False
            self.method3Collapsible.visible = False
            # Uncheck other buttons
            self.roiMethodButton.setChecked(False)
            self.brushMethodButton.setChecked(False)
        else:
            # Hide Click & Grow method panel
            self.method2Collapsible.visible = False

    def onBrushMethodSelected(self, checked):
        """Toggle Brush method panel"""
        if checked:
            # Show Brush method panel, hide others
            self.method1Collapsible.visible = False
            self.method2Collapsible.visible = False
            self.method3Collapsible.visible = True
            self.method3Collapsible.collapsed = False
            # Uncheck other buttons
            self.roiMethodButton.setChecked(False)
            self.clickGrowMethodButton.setChecked(False)
        else:
            # Hide Brush method panel
            self.method3Collapsible.visible = False

    # ============================================================
    # Segmentation Callbacks
    # ============================================================

    def onPlaceROIToggled(self, checked):
        """Toggle ROI placement"""
        if checked:
            if not self.roiNode:
                self.roiNode = self.logic.createROI(self.volumeNode)
                print(f"ROI created: {self.roiNode.GetName()} at position {self.roiNode.GetXYZ([0,0,0])}")
            self.logic.startEditingROI(self.roiNode)
            self.logic.setROIVisibility(self.roiNode, True)
            self.placeROIButton.text = "Finish ROI"
            self.applyROIThresholdButton.enabled = True
            self.showROICheckbox.enabled = True
            self.showROICheckbox.setChecked(True)

            # Enable rotation controls
            self.roiRotationLR.enabled = True
            self.roiRotationPA.enabled = True
            self.roiRotationIS.enabled = True
            self.resetRotationButton.enabled = True

            # Switch to 3D view to see the ROI
            layoutManager = slicer.app.layoutManager()
            threeDWidget = layoutManager.threeDWidget(0)
            threeDView = threeDWidget.threeDView()
            threeDView.resetFocalPoint()

            slicer.util.infoDisplay(
                "ROI Box placed!\n\n"
                "You can now:\n"
                "1. Drag the box to position it over the aortic valve\n"
                "2. Resize using the handles\n"
                "3. Use the rotation sliders (L-R, P-A, I-S) to rotate the box\n"
                "4. Click 'Finish ROI' when ready\n"
                "5. Then click 'Apply Threshold in ROI'"
            )
        else:
            self.logic.stopEditingROI()
            self.placeROIButton.text = "Place ROI Box"

    def onApplyROIThreshold(self):
        """Apply threshold within ROI and auto-hide ROI"""
        if not self.segmentationNode or not self.roiNode:
            slicer.util.errorDisplay("Please place an ROI first")
            return

        # Use standard non-contrast threshold (130 HU)
        thresholdToUse = self.thresholdValue

        self.logic.applyThresholdInROI(
            self.volumeNode,
            self.segmentationNode,
            self.roiNode,
            thresholdToUse
        )

        # Auto-hide ROI after applying threshold
        self.logic.setROIVisibility(self.roiNode, False)
        self.showROICheckbox.setChecked(False)

        # Enable calculation button and Results tab
        self.calculateButton.enabled = True
        self.tabWidget.setTabEnabled(3, True)  # Enable Results tab

        slicer.util.infoDisplay(
            f"Segmentation applied with threshold {thresholdToUse} HU (Non-Contrast CT)\n\n"
            "ROI has been hidden. Use 'Show ROI' checkbox to display it again if needed.\n\n"
            "You can now proceed to Tab 3 (Results & Analysis) to calculate the Agatston Score."
        )

    def onToggleROIVisibility(self, checked):
        """Toggle ROI visibility"""
        if self.roiNode:
            self.logic.setROIVisibility(self.roiNode, checked)

    def onROIRotationChanged(self, value):
        """Handle ROI rotation slider changes"""
        if not self.roiNode:
            return

        # Get rotation values from sliders
        rotLR = self.roiRotationLR.value
        rotPA = self.roiRotationPA.value
        rotIS = self.roiRotationIS.value

        # Apply rotation using transform
        self.logic.setROIRotation(self.roiNode, rotLR, rotPA, rotIS)

    def onResetROIRotation(self):
        """Reset ROI rotation to zero"""
        if not self.roiNode:
            return

        # Reset sliders
        self.roiRotationLR.value = 0
        self.roiRotationPA.value = 0
        self.roiRotationIS.value = 0

        # Reset transform
        self.logic.setROIRotation(self.roiNode, 0, 0, 0)

    # ============================================================================
    # CLICK & GROW MODE CALLBACKS
    # ============================================================================

    def onClickGrowToggled(self, checked):
        """Toggle click and grow mode"""
        if checked:
            # Disable other modes
            self.paintButton.setChecked(False)
            self.eraseButton.setChecked(False)

            # Start click and grow mode
            self.logic.startClickGrowMode(self.segmentationNode, self.volumeNode, self.thresholdValue)
            self.clickGrowButton.text = "Stop Click & Grow"

            # Enable calculation button and Results tab
            self.calculateButton.enabled = True
            self.tabWidget.setTabEnabled(3, True)  # Enable Results tab

            slicer.util.infoDisplay(
                "Click & Grow Mode enabled!\n\n"
                "Click on any calcification in the 2D slices or 3D view.\n"
                "The algorithm will automatically segment the calcification\n"
                "and propagate to adjacent slices in 3D.\n\n"
                "Only voxels ≥ 130 HU will be included.\n\n"
                "Click 'Stop Click & Grow' when finished."
            )
        else:
            self.logic.stopClickGrowMode()
            self.clickGrowButton.text = "Enable Click & Grow"

    # ============================================================================
    # PAINT/ERASE MODE CALLBACKS
    # ============================================================================

    def onPaintToggled(self, checked):
        """Toggle paint mode"""
        if checked:
            # Disable other modes
            self.clickGrowButton.setChecked(False)
            self.eraseButton.setChecked(False)

            # Use standard non-contrast threshold (130 HU)
            thresholdToUse = self.thresholdValue

            brushSize = self.brushSizeSlider.value
            self.logic.startPaintMode(self.segmentationNode, thresholdToUse, brushSize)
            self.paintButton.text = "Stop Painting"
            # Enable calculation button and Results tab after painting starts
            self.calculateButton.enabled = True
            self.tabWidget.setTabEnabled(3, True)  # Enable Results tab
        else:
            self.logic.stopPaintMode()
            self.paintButton.text = "Paint Mode"

    def onEraseToggled(self, checked):
        """Toggle erase mode"""
        if checked:
            # Disable other modes
            self.clickGrowButton.setChecked(False)
            self.paintButton.setChecked(False)

            brushSize = self.brushSizeSlider.value
            self.logic.startEraseMode(self.segmentationNode, brushSize)
            self.eraseButton.text = "Stop Erasing"
        else:
            self.logic.stopEraseMode()
            self.eraseButton.text = "Erase Mode"

    def onClearSegmentation(self):
        """Clear segmentation"""
        if self.segmentationNode:
            self.logic.clearSegmentation(self.segmentationNode)
            self.calculateButton.enabled = False

    def onPreview3DToggled(self, checked):
        """Toggle 3D preview"""
        if self.segmentationNode:
            displayNode = self.segmentationNode.GetDisplayNode()
            if displayNode:
                displayNode.SetVisibility(checked)

    def onCalculate(self):
        """Calculate Agatston score and enable report tab"""
        if not self.segmentationNode or not self.volumeNode:
            slicer.util.errorDisplay("No segmentation available")
            return

        # Get patient sex and age for classification
        sex = 'F' if self.femaleButton.isChecked() else 'M'
        age = None
        try:
            ageStr = self.patientAgeEdit.text
            if ageStr:
                age = int(ageStr)
        except:
            pass

        # Calculate
        results = self.logic.calculateAgatstonScore(
            self.volumeNode,
            self.segmentationNode,
            sex,
            age
        )

        # Update UI
        self.updateResultsDisplay(results)

        # Enable visualization buttons and report tab
        self.show3DButton.enabled = True
        self.showChartsButton.enabled = True
        self.showBullsEyeButton.enabled = True
        self.tabWidget.setTabEnabled(4, True)  # Enable Report tab (index 4 with Settings tab)
        self.generateReportButton.enabled = True
        self.reportStatusLabel.setText("Status: Ready to generate report!")
        self.reportStatusLabel.setStyleSheet("color: #28a745; font-weight: bold; padding: 5px;")

        slicer.util.infoDisplay(
            "Calculation complete!\n\n"
            f"Agatston Score: {results['agatston_score']:.1f} AU\n"
            f"Classification: {results['classification']}\n\n"
            "View detailed results below and generate report in Tab 4 (Generate Report)."
        )

    def updateResultsDisplay(self, results):
        """Update results display"""
        # Main score
        agatston = results['agatston_score']
        self.agatstonScoreLabel.text = f"Agatston Score: {agatston:.1f} AU\n(Non-Contrast CT - Standard Agatston)"

        # Classification
        classification = results['classification']
        severity = results['severity']

        colorMap = {
            'Normal/Minimal': '#4CAF50',
            'Mild': '#FFC107',
            'Moderate': '#FF9800',
            'Severe': '#F44336'
        }

        color = colorMap.get(severity, '#cccccc')

        # Add severity threshold info
        sex = results.get('patient_sex', 'M')
        severeThreshold = 1300 if sex == 'F' else 2000
        thresholdInfo = f"<br><small>Non-contrast CT threshold: ≥{severeThreshold} AU for severe AS</small>"

        self.classificationLabel.setText(f"Classification: {classification}{thresholdInfo}")
        self.classificationLabel.setStyleSheet(f"""
            QLabel {{
                font-size: 16px;
                font-weight: bold;
                padding: 10px;
                background-color: {color};
                color: white;
                border-radius: 5px;
            }}
        """)

        # Table
        metrics = [
            ("Total Volume (mm³)", f"{results['total_volume_mm3']:.2f}"),
            ("Equivalent Mass (mg)", f"{results['equivalent_mass_mg']:.2f}"),
            ("Number of Lesions", str(results['num_lesions'])),
            ("Mean Density (HU)", f"{results['mean_density']:.1f}"),
            ("Max Density (HU)", f"{results['max_density']:.1f}"),
            ("Density Factor 1 (130-199 HU)", f"{results['density_distribution'][0]:.1f}%"),
            ("Density Factor 2 (200-299 HU)", f"{results['density_distribution'][1]:.1f}%"),
            ("Density Factor 3 (300-399 HU)", f"{results['density_distribution'][2]:.1f}%"),
            ("Density Factor 4 (≥400 HU)", f"{results['density_distribution'][3]:.1f}%"),
        ]

        self.resultsTable.setRowCount(len(metrics))
        for i, (metric, value) in enumerate(metrics):
            self.resultsTable.setItem(i, 0, qt.QTableWidgetItem(metric))
            self.resultsTable.setItem(i, 1, qt.QTableWidgetItem(value))

    def onShow3D(self):
        """Show 3D model"""
        if self.segmentationNode:
            colorByDensity = self.colorByDensityCheckbox.isChecked()
            self.logic.create3DModel(self.segmentationNode, self.volumeNode, colorByDensity)

    def onShowCharts(self):
        """Show charts"""
        if self.segmentationNode and self.volumeNode:
            self.logic.createCharts(self.volumeNode, self.segmentationNode)

    def onShowBullsEye(self):
        """Show bull's eye plot"""
        if self.segmentationNode and self.volumeNode:
            self.logic.createBullsEyeAnalysis(self.volumeNode, self.segmentationNode)

    def onGenerateReport(self):
        """Generate PDF report"""
        outputPath = self.outputPathEdit.text
        if not outputPath:
            outputPath = qt.QFileDialog.getExistingDirectory(None, "Select Output Directory")
            if not outputPath:
                return
            self.outputPathEdit.text = outputPath

        # Update patient info
        self.updatePatientInfo()

        # Generate report
        reportOptions = {
            'includeScreenshots': self.includeScreenshotsCheckbox.isChecked(),
            'include3D': self.include3DCheckbox.isChecked(),
            'includeCharts': self.includeChartsCheckbox.isChecked()
        }

        success = self.logic.generatePDFReport(
            self.volumeNode,
            self.segmentationNode,
            self.patientInfo,
            outputPath,
            reportOptions
        )

        if success:
            slicer.util.infoDisplay(f"Report generated successfully!\nSaved to: {outputPath}")
        else:
            slicer.util.errorDisplay("Failed to generate report")

    def onBrowseOutput(self):
        """Browse for output directory"""
        directory = qt.QFileDialog.getExistingDirectory(None, "Select Output Directory")
        if directory:
            self.outputPathEdit.text = directory

    def updatePatientInfo(self):
        """Update patient information dictionary"""
        self.patientInfo['name'] = self.patientNameEdit.text
        self.patientInfo['id'] = self.patientIDEdit.text
        self.patientInfo['sex'] = 'F' if self.femaleButton.isChecked() else 'M'
        self.patientInfo['age'] = self.patientAgeEdit.text
        self.patientInfo['date'] = datetime.now().strftime("%Y-%m-%d %H:%M")

#
# SlicerCaScoreLogic
#

class AoCaScoreLogic(ScriptedLoadableModuleLogic):
    """Logic for aortic valve calcium scoring"""

    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self.currentResults = {}

    def setupOptimalView(self, volumeNode, applyWindowLevel=True):
        """Setup optimal view layout for aortic valve"""
        # Set layout to Four-Up
        layoutManager = slicer.app.layoutManager()
        layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)

        # Configure each slice view
        sliceNames = ['Red', 'Yellow', 'Green']
        orientations = ['Axial', 'Sagittal', 'Coronal']

        for sliceName, orientation in zip(sliceNames, orientations):
            sliceWidget = layoutManager.sliceWidget(sliceName)
            sliceNode = sliceWidget.sliceLogic().GetSliceNode()

            # Set orientation
            sliceNode.SetOrientation(orientation)

            # Set background volume
            sliceWidget.sliceLogic().GetSliceCompositeNode().SetBackgroundVolumeID(volumeNode.GetID())

            if applyWindowLevel:
                # Apply calcium window/level preset
                displayNode = volumeNode.GetDisplayNode()
                if displayNode:
                    displayNode.SetWindow(1500)
                    displayNode.SetLevel(300)

            # Reset field of view
            sliceWidget.sliceLogic().FitSliceToAll()

        # Setup 3D view
        threeDWidget = layoutManager.threeDWidget(0)
        threeDView = threeDWidget.threeDView()
        threeDView.resetFocalPoint()

    def enableSliceIntersections(self, enable=True):
        """Enable or disable slice intersection visibility"""
        # Get all slice nodes
        sliceNodes = slicer.util.getNodesByClass('vtkMRMLSliceNode')

        for sliceNode in sliceNodes:
            # Get or create slice display node
            displayNode = sliceNode.GetSliceDisplayNode()
            if not displayNode:
                sliceNode.CreateDefaultDisplayNodes()
                displayNode = sliceNode.GetSliceDisplayNode()

            if displayNode:
                if enable:
                    # Enable intersecting slices visibility
                    displayNode.SetIntersectingSlicesVisibility(True)
                    displayNode.SetIntersectingSlicesInteractive(True)
                else:
                    displayNode.SetIntersectingSlicesVisibility(False)
                    displayNode.SetIntersectingSlicesInteractive(False)

        # Force refresh of all views
        slicer.app.processEvents()

    def createSegmentation(self, volumeNode):
        """Create new segmentation for calcium"""
        # Remove any existing segmentation with the same name to avoid conflicts
        existingNodes = slicer.util.getNodesByClass("vtkMRMLSegmentationNode")
        for node in existingNodes:
            if node.GetName() == "AorticValveCalcium":
                print(f"Removing old segmentation node: {node.GetName()}")
                slicer.mrmlScene.RemoveNode(node)

        # Create segmentation node and ensure it's in the scene
        segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "AorticValveCalcium")

        # Set reference geometry BEFORE creating display nodes
        segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)

        # Create calcium segment first
        segmentation = segmentationNode.GetSegmentation()
        segmentId = segmentation.AddEmptySegment("Calcium")
        segment = segmentation.GetSegment(segmentId)

        # Set color
        segment.SetColor(1.0, 0.0, 0.0)  # Red

        # Create display nodes
        if segmentationNode.GetScene():
            segmentationNode.CreateDefaultDisplayNodes()
        else:
            print("Warning: Segmentation node not in scene, cannot create display nodes")

        # Ensure subject hierarchy is updated
        shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
        if shNode:
            itemID = shNode.GetItemByDataNode(segmentationNode)
            if not itemID:
                itemID = shNode.CreateItem(shNode.GetSceneItemID(), segmentationNode)
                print(f"✓ Created subject hierarchy item for segmentation: {itemID}")

        print(f"✓ Segmentation node created: {segmentationNode.GetName()}")
        return segmentationNode

    def createROI(self, volumeNode):
        """Create ROI for valve region"""
        roiNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode", "ValveROI")
        roiNode.CreateDefaultDisplayNodes()

        # Set initial size based on volume
        bounds = [0] * 6
        volumeNode.GetBounds(bounds)
        center = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]

        # Default ROI size for aortic valve (approximately 40x40x40 mm)
        roiNode.SetXYZ(center)
        roiNode.SetRadiusXYZ(20, 20, 20)

        return roiNode

    def startEditingROI(self, roiNode):
        """Start editing ROI"""
        selectionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLSelectionNodeSingleton")
        selectionNode.SetActivePlaceNodeID(roiNode.GetID())
        interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
        interactionNode.SetCurrentInteractionMode(interactionNode.Place)

    def stopEditingROI(self):
        """Stop editing ROI"""
        interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
        interactionNode.SetCurrentInteractionMode(interactionNode.ViewTransform)

    def setROIVisibility(self, roiNode, visible):
        """Set ROI visibility"""
        if roiNode:
            displayNode = roiNode.GetDisplayNode()
            if displayNode:
                displayNode.SetVisibility(visible)

    def setROIRotation(self, roiNode, rotLR, rotPA, rotIS):
        """Set ROI rotation using transform

        Args:
            roiNode: The ROI node to rotate
            rotLR: Rotation around Left-Right axis (degrees)
            rotPA: Rotation around Posterior-Anterior axis (degrees)
            rotIS: Rotation around Inferior-Superior axis (degrees)
        """
        if not roiNode:
            return

        # Get or create transform node
        transformNode = roiNode.GetParentTransformNode()
        if not transformNode:
            transformNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode", "ROI_Transform")
            roiNode.SetAndObserveTransformNodeID(transformNode.GetID())

        # Create rotation matrix
        import vtk
        transform = vtk.vtkTransform()

        # Get ROI center to rotate around it
        center = [0, 0, 0]
        roiNode.GetXYZ(center)

        # Apply rotations in RAS coordinate system
        # Move to origin, rotate, move back
        transform.Translate(center[0], center[1], center[2])
        transform.RotateX(rotLR)  # L-R axis
        transform.RotateY(rotPA)  # P-A axis
        transform.RotateZ(rotIS)  # I-S axis
        transform.Translate(-center[0], -center[1], -center[2])

        # Apply transform to node
        transformNode.SetMatrixTransformToParent(transform.GetMatrix())

    # ============================================================================
    # ROI-BASED SEGMENTATION METHODS (Valve ROI only)
    # ============================================================================

    def applyThresholdInROI(self, volumeNode, segmentationNode, roiNode, threshold):
        """Apply threshold segmentation within ROI"""
        # Step 1: Crop volume with ROI to create a temporary volume
        cropVolumeLogic = slicer.modules.cropvolume.logic()
        cropVolumeParameterNode = slicer.vtkMRMLCropVolumeParametersNode()
        cropVolumeParameterNode.SetROINodeID(roiNode.GetID())
        cropVolumeParameterNode.SetInputVolumeNodeID(volumeNode.GetID())
        cropVolumeParameterNode.SetVoxelBased(True)

        slicer.mrmlScene.AddNode(cropVolumeParameterNode)

        # Create output volume
        croppedVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "TempCroppedVolume")
        cropVolumeParameterNode.SetOutputVolumeNodeID(croppedVolume.GetID())

        # Apply crop
        cropVolumeLogic.Apply(cropVolumeParameterNode)

        # Step 2: Apply threshold on the cropped volume
        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
        segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
        segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
        segmentEditorWidget.setSegmentationNode(segmentationNode)

        # Temporarily set the cropped volume as source
        segmentEditorWidget.setSourceVolumeNode(croppedVolume)

        # Get the calcium segment
        calciumSegmentId = segmentationNode.GetSegmentation().GetNthSegmentID(0)
        segmentEditorWidget.setCurrentSegmentID(calciumSegmentId)

        # Apply threshold
        segmentEditorWidget.setActiveEffectByName("Threshold")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("MinimumThreshold", str(threshold))
        effect.setParameter("MaximumThreshold", "3000")
        effect.self().onApply()

        # Step 3: Restore original volume as source
        segmentEditorWidget.setSourceVolumeNode(volumeNode)

        # Step 4: Clean up temporary nodes
        slicer.mrmlScene.RemoveNode(croppedVolume)
        slicer.mrmlScene.RemoveNode(cropVolumeParameterNode)
        slicer.mrmlScene.RemoveNode(segmentEditorNode)

    def startPaintMode(self, segmentationNode, threshold, brushSize):
        """Start threshold paint mode"""
        # Switch to Segment Editor
        slicer.util.selectModule("SegmentEditor")

        # Get segment editor widget
        segmentEditorWidget = slicer.modules.segmenteditor.widgetRepresentation().self().editor
        segmentEditorWidget.setSegmentationNode(segmentationNode)

        # Select segment
        segmentId = segmentationNode.GetSegmentation().GetNthSegmentID(0)
        segmentEditorWidget.setCurrentSegmentID(segmentId)

        # Activate Paint effect
        segmentEditorWidget.setActiveEffectByName("Paint")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("BrushRelativeDiameter", str(brushSize))

        # Set masking to threshold
        segmentEditorNode = segmentEditorWidget.mrmlSegmentEditorNode()
        segmentEditorNode.SetMasterVolumeIntensityMask(True)
        segmentEditorNode.SetMasterVolumeIntensityMaskRange(threshold, 3000)

    def stopPaintMode(self):
        """Stop paint mode"""
        # Return to this module
        slicer.util.selectModule("SlicerCaScore")

    def startEraseMode(self, segmentationNode, brushSize):
        """Start erase mode"""
        slicer.util.selectModule("SegmentEditor")

        segmentEditorWidget = slicer.modules.segmenteditor.widgetRepresentation().self().editor
        segmentEditorWidget.setSegmentationNode(segmentationNode)

        segmentId = segmentationNode.GetSegmentation().GetNthSegmentID(0)
        segmentEditorWidget.setCurrentSegmentID(segmentId)

        segmentEditorWidget.setActiveEffectByName("Erase")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("BrushRelativeDiameter", str(brushSize))

    def stopEraseMode(self):
        """Stop erase mode"""
        slicer.util.selectModule("SlicerCaScore")

    def startClickGrowMode(self, segmentationNode, volumeNode, threshold):
        """Start click and grow mode with automatic 3D region growing on click"""
        self.clickGrowSegmentationNode = segmentationNode
        self.clickGrowVolumeNode = volumeNode
        self.clickGrowThreshold = threshold
        self.clickGrowObserverTag = None

        # Store volume array for fast access
        self.clickGrowVolumeArray = slicer.util.arrayFromVolume(volumeNode)

        # Set up markup fiducial for clicking
        self.clickGrowFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", "ClickGrow_Seeds")
        self.clickGrowFiducialNode.CreateDefaultDisplayNodes()

        # Make fiducials small and temporary
        displayNode = self.clickGrowFiducialNode.GetDisplayNode()
        displayNode.SetGlyphScale(1.5)
        displayNode.SetSelectedColor(0, 1, 0)  # Green

        # Set interaction mode to place fiducials
        selectionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLSelectionNodeSingleton")
        selectionNode.SetActivePlaceNodeID(self.clickGrowFiducialNode.GetID())
        interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
        interactionNode.SetCurrentInteractionMode(interactionNode.Place)

        # Add observer to detect when fiducial is placed
        self.clickGrowObserverTag = self.clickGrowFiducialNode.AddObserver(
            slicer.vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self.onFiducialPlaced
        )

        print(f"Click & Grow Mode: Active with threshold ≥ {threshold} HU")
        print("Click on any calcification to segment it automatically")

    def onFiducialPlaced(self, caller, event):
        """Called when user clicks to place a fiducial - perform region growing"""
        try:
            fiducialNode = caller
            numFiducials = fiducialNode.GetNumberOfControlPoints()

            if numFiducials == 0:
                return

            # Get the position of the last placed fiducial
            lastIndex = numFiducials - 1
            seedPos_RAS = [0, 0, 0]
            fiducialNode.GetNthControlPointPosition(lastIndex, seedPos_RAS)

            print(f"\nClick #{lastIndex + 1}: Growing from seed at RAS: {seedPos_RAS}")

            # Convert RAS to IJK
            volumeNode = self.clickGrowVolumeNode
            rasToIjkMatrix = vtk.vtkMatrix4x4()
            volumeNode.GetRASToIJKMatrix(rasToIjkMatrix)

            seedPos_IJK = [0, 0, 0, 1]
            rasPoint = [seedPos_RAS[0], seedPos_RAS[1], seedPos_RAS[2], 1]
            ijkPoint = rasToIjkMatrix.MultiplyPoint(rasPoint)
            seedPos_IJK = [int(round(ijkPoint[0])), int(round(ijkPoint[1])), int(round(ijkPoint[2]))]

            print(f"Seed IJK: {seedPos_IJK}")

            # Perform 3D region growing from this seed
            self.performRegionGrowing3D(seedPos_IJK)

            # Simply remove the fiducial point after processing
            # Keep the same node to avoid complexity
            fiducialNode.RemoveNthControlPoint(lastIndex)

            # Force placement mode to stay active
            qt.QTimer.singleShot(100, self.reactivatePlacementMode)

        except Exception as e:
            print(f"Error in onFiducialPlaced: {str(e)}")
            import traceback
            traceback.print_exc()

    def reactivatePlacementMode(self):
        """Reactivate placement mode for continuous clicking"""
        try:
            if hasattr(self, 'clickGrowFiducialNode') and self.clickGrowFiducialNode:
                selectionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLSelectionNodeSingleton")
                selectionNode.SetActivePlaceNodeID(self.clickGrowFiducialNode.GetID())
                interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
                interactionNode.SetCurrentInteractionMode(interactionNode.Place)
                print("→ Ready for next click...")
        except Exception as e:
            print(f"Error reactivating placement mode: {str(e)}")

    def performRegionGrowing3D(self, seedIJK):
        """Perform 3D region growing from seed point using simple connected components"""
        volumeArray = self.clickGrowVolumeArray
        threshold = self.clickGrowThreshold

        # Check if seed is valid
        dims = volumeArray.shape
        i, j, k = seedIJK

        if not (0 <= k < dims[0] and 0 <= j < dims[1] and 0 <= i < dims[2]):
            print(f"Seed outside volume bounds: {seedIJK}")
            return

        seedValue = volumeArray[k, j, i]
        print(f"Seed HU value: {seedValue:.1f}")

        if seedValue < threshold:
            print(f"Seed HU ({seedValue:.1f}) is below threshold ({threshold}). Not segmenting.")
            slicer.util.warningDisplay(
                f"Clicked point has HU = {seedValue:.1f}\n"
                f"This is below the threshold ({threshold} HU).\n\n"
                f"Please click on a calcification (bright white area)."
            )
            return

        # Simple and effective: 3D connected component labeling with threshold
        from scipy import ndimage

        print(f"Performing 3D region growing from seed [{i}, {j}, {k}]...")

        # Create binary mask: voxels >= threshold
        binaryMask = (volumeArray >= threshold).astype(np.uint8)

        # Label all connected components in 3D (26-connectivity for full 3D)
        # This finds all calcifications in the entire volume
        labeled, numFeatures = ndimage.label(binaryMask, structure=ndimage.generate_binary_structure(3, 3))

        # Find which label the seed belongs to
        seedLabel = labeled[k, j, i]

        if seedLabel == 0:
            print("Seed is not part of any connected component above threshold")
            slicer.util.warningDisplay(
                f"No connected calcification found at this point.\n\n"
                f"The clicked voxel may be isolated or below threshold."
            )
            return

        # Extract only the connected component containing the seed
        # This is the "click and grow" - we only take the lesion we clicked on
        componentMask = (labeled == seedLabel).astype(np.uint8)

        # Count voxels
        numVoxels = np.sum(componentMask)

        # Calculate volume
        volumeNode = self.clickGrowVolumeNode
        spacing = volumeNode.GetSpacing()
        voxelVolume = spacing[0] * spacing[1] * spacing[2]  # mm³
        lesionVolume = numVoxels * voxelVolume

        print(f"✓ Segmented lesion: {numVoxels} voxels, {lesionVolume:.1f} mm³")

        # Add to main segmentation
        segmentationNode = self.clickGrowSegmentationNode
        segmentation = segmentationNode.GetSegmentation()
        mainSegmentId = segmentation.GetNthSegmentID(0)

        # Ensure the segmentation has the correct geometry
        # This is crucial for proper export/import
        segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)

        # Check if segment has any data yet
        # For the first click, the segment is empty and export will fail
        segment = segmentation.GetSegment(mainSegmentId)

        # Try to get current segmentation array
        currentArray = None
        if segment and segment.GetRepresentation("Binary labelmap"):
            # Segment has data, try to export it
            try:
                currentArray = slicer.util.arrayFromSegmentBinaryLabelmap(
                    segmentationNode, mainSegmentId, volumeNode
                )
            except Exception as e:
                print(f"Note: Could not export segment (first click?): {e}")

        # If no current array (first click or export failed), create empty one
        if currentArray is None:
            imageData = volumeNode.GetImageData()
            dims = imageData.GetDimensions()
            currentArray = np.zeros((dims[2], dims[1], dims[0]), dtype=np.uint8)
            print("✓ Initializing segmentation with first calcification")

        # Add new component (OR operation) - accumulates multiple clicks
        currentArray = np.logical_or(currentArray, componentMask).astype(np.uint8)

        # Update main segmentation
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            currentArray,
            segmentationNode,
            mainSegmentId,
            volumeNode
        )

        print(f"✓ Successfully added calcification to segmentation!")

        # Show brief notification
        slicer.util.showStatusMessage(f"Added lesion: {lesionVolume:.1f} mm³", 2000)

    def stopClickGrowMode(self):
        """Stop click and grow mode"""
        # Remove observer
        if hasattr(self, 'clickGrowFiducialNode') and self.clickGrowFiducialNode and self.clickGrowObserverTag:
            self.clickGrowFiducialNode.RemoveObserver(self.clickGrowObserverTag)
            self.clickGrowObserverTag = None

        # Remove fiducial node
        if hasattr(self, 'clickGrowFiducialNode') and self.clickGrowFiducialNode:
            slicer.mrmlScene.RemoveNode(self.clickGrowFiducialNode)
            self.clickGrowFiducialNode = None

        # Stop placing mode
        interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
        interactionNode.SetCurrentInteractionMode(interactionNode.ViewTransform)

        # Clean up references
        self.clickGrowSegmentationNode = None
        self.clickGrowVolumeNode = None
        self.clickGrowThreshold = None
        self.clickGrowVolumeArray = None

        print("Click & Grow Mode: Stopped")

    def clearSegmentation(self, segmentationNode):
        """Clear all segments"""
        segmentation = segmentationNode.GetSegmentation()
        segmentId = segmentation.GetNthSegmentID(0)

        if segmentId:
            slicer.vtkSlicerSegmentationsModuleLogic.ClearSegment(segmentationNode, segmentId)

    def calculateAgatstonScore(self, volumeNode, segmentationNode, patientSex, patientAge=None):
        """Calculate Agatston score and related metrics (Non-Contrast CT only)

        Args:
            volumeNode: CT volume node
            segmentationNode: Segmentation node containing calcium
            patientSex: 'M' or 'F' for sex-specific classification
            patientAge: Patient age (optional, for percentile comparison)
        """
        # Get segment as binary labelmap
        segmentId = segmentationNode.GetSegmentation().GetNthSegmentID(0)
        segmentArray = slicer.util.arrayFromSegmentBinaryLabelmap(segmentationNode, segmentId, volumeNode)
        volumeArray = slicer.util.arrayFromVolume(volumeNode)

        if segmentArray is None or not np.any(segmentArray):
            self.currentResults = self.getEmptyResults()
            return self.currentResults

        # Get spacing for area/volume calculation
        spacing = volumeNode.GetSpacing()
        sliceArea = spacing[0] * spacing[1]  # mm² per pixel in slice
        voxelVolume = spacing[0] * spacing[1] * spacing[2]  # mm³ per voxel

        # Try to get slice thickness from DICOM metadata
        sliceThickness = None
        try:
            # Try to read SliceThickness from DICOM
            import pydicom
            # Get the file path from the volume node
            storageNode = volumeNode.GetStorageNode()
            if storageNode:
                filePath = storageNode.GetFileName()
                if filePath:
                    dicom = pydicom.dcmread(filePath, stop_before_pixels=True)
                    if hasattr(dicom, 'SliceThickness'):
                        sliceThickness = float(dicom.SliceThickness)
                        print(f"ℹ DICOM Slice Thickness: {sliceThickness:.1f}mm")
        except:
            # If we can't read DICOM metadata, assume thickness = spacing
            pass

        # If we couldn't read thickness from DICOM, assume it equals spacing
        if sliceThickness is None:
            sliceThickness = spacing[2]
            print(f"ℹ Slice thickness not found in DICOM, assuming thickness = spacing = {sliceThickness:.1f}mm")

        # Label connected components (lesions)
        # IMPORTANT: Original Agatston uses 2D connectivity (per-slice labeling)
        # This matches commercial software like Syngo.via
        from scipy import ndimage

        # Standard Agatston uses 3.0mm slice thickness
        STANDARD_SLICE_THICKNESS = 3.0  # mm (Agatston standard)
        sliceSpacing = spacing[2]

        # Normalization factor for non-3mm slices
        # Process ALL slices and normalize score at the end
        # Reference: Lassoan et al., 3D Slicer Agatston implementation
        sliceNormalizationFactor = sliceSpacing / STANDARD_SLICE_THICKNESS
        sliceStep = 1  # Process all slices

        print(f"Slice spacing: {sliceSpacing:.2f}mm, normalization factor: {sliceNormalizationFactor:.3f}")
        if abs(sliceSpacing - STANDARD_SLICE_THICKNESS) < 0.1:
            print(f"  -> Standard 3mm acquisition, no normalization needed")
        else:
            print(f"  -> Score will be normalized to 3mm equivalent")

        # Use 2D slice-by-slice labeling (Agatston standard)
        # This prevents overestimation from 3D connectivity
        labeled_array = np.zeros_like(segmentArray, dtype=np.int32)
        lesion_counter = 0

        for sliceIdx in range(0, segmentArray.shape[0], sliceStep):  # Skip slices if overlapping
            sliceMask = segmentArray[sliceIdx, :, :]
            if np.any(sliceMask):
                # Label lesions in this slice only (2D connectivity)
                labeled_slice, n_lesions_in_slice = ndimage.label(sliceMask)
                # Offset labels to make them unique across all slices
                labeled_slice[labeled_slice > 0] += lesion_counter
                labeled_array[sliceIdx, :, :] = labeled_slice
                lesion_counter += n_lesions_in_slice

        num_lesions = lesion_counter
        print(f"2D slice-by-slice labeling: {num_lesions} lesions detected (step={sliceStep})")

        # Agatston standard: minimum requirements PER SLICE
        # Clinical standard: minimum 1 mm2 area per slice (Agatston et al. 1990)
        # Reference: "Only contiguous voxels totaling >=1 mm2 are counted as lesions"
        MIN_AREA_MM2_STANDARD = 1.0  # mm2 - fixed clinical standard

        # Calculate minimum pixels needed to reach 1 mm2 threshold
        # Use at least 2 pixels to avoid single-pixel noise, but ensure >=1 mm2
        minPixelsFor1mm2 = int(np.ceil(MIN_AREA_MM2_STANDARD / sliceArea))
        MIN_PIXELS_PER_SLICE = max(2, minPixelsFor1mm2)  # At least 2 pixels, but >=1 mm2
        MIN_AREA_MM2 = MIN_PIXELS_PER_SLICE * sliceArea  # Actual threshold applied

        # Print spacing information for debugging
        print(f"Image spacing: {spacing[0]:.3f} x {spacing[1]:.3f} x {spacing[2]:.3f} mm")
        print(f"Slice area per pixel: {sliceArea:.3f} mm2")
        print(f"Voxel volume: {voxelVolume:.3f} mm3")
        print(f"Minimum pixels required for 1 mm2 standard: {minPixelsFor1mm2}")
        print(f"Minimum area threshold per slice: {MIN_AREA_MM2:.3f} mm2 ({MIN_PIXELS_PER_SLICE} pixels)")

        # Verify compliance with 1 mm2 standard
        if MIN_AREA_MM2 < MIN_AREA_MM2_STANDARD:
            print(f"WARNING: Threshold {MIN_AREA_MM2:.3f} mm2 is below clinical standard of {MIN_AREA_MM2_STANDARD} mm2")
        else:
            print(f"OK: Threshold complies with Agatston standard (>={MIN_AREA_MM2_STANDARD} mm2)")

        # Calculate metrics
        totalAgatston = 0
        totalVolume = 0
        densityDistribution = [0, 0, 0, 0]  # Count for each density factor
        allDensities = []
        validLesionCount = 0
        totalSlicesProcessed = 0
        totalSlicesFiltered = 0

        for lesion_id in range(1, num_lesions + 1):
            lesionMask = (labeled_array == lesion_id)
            lesionVoxels = volumeArray[lesionMask]

            # Get all densities for statistics
            allDensities.extend(lesionVoxels.tolist())

            # Calculate Agatston score slice-by-slice (STANDARD METHOD)
            # Density factor must be determined PER SLICE, not per lesion!
            lesionScore = 0
            lesionVolume = 0
            lesionHasValidSlice = False

            for sliceIdx in range(lesionMask.shape[0]):
                sliceMask = lesionMask[sliceIdx, :, :]
                if np.any(sliceMask):
                    # Calculate area for this slice
                    pixelCount = np.sum(sliceMask)
                    area_mm2 = pixelCount * sliceArea

                    # CRITICAL: Filter by minimum area PER SLICE (Agatston standard)
                    if area_mm2 < MIN_AREA_MM2:
                        totalSlicesFiltered += 1
                        continue  # Skip this slice if too small

                    totalSlicesProcessed += 1
                    lesionHasValidSlice = True

                    # Get max density in THIS SLICE of the lesion
                    sliceVoxels = volumeArray[sliceIdx, :, :][sliceMask]
                    maxDensityInSlice = np.max(sliceVoxels)

                    # Determine density factor for THIS SLICE (Agatston standard)
                    if maxDensityInSlice >= 400:
                        densityFactor = 4
                    elif maxDensityInSlice >= 300:
                        densityFactor = 3
                    elif maxDensityInSlice >= 200:
                        densityFactor = 2
                    else:  # 130-199
                        densityFactor = 1

                    # Calculate score for this slice
                    lesionScore += area_mm2 * densityFactor

                    # Volume calculation must account for slice skipping
                    # If we're using every 2nd slice, each slice represents sliceStep slices
                    # Volume per slice = pixelCount * sliceArea * (sliceSpacing * sliceStep)
                    sliceVolumeContribution = pixelCount * sliceArea * (sliceSpacing * sliceStep)
                    lesionVolume += sliceVolumeContribution

            # Only count lesion if it has at least one valid slice
            if lesionHasValidSlice:
                validLesionCount += 1
                totalAgatston += lesionScore
                totalVolume += lesionVolume

                # Track density distribution (use max density of entire lesion for statistics)
                maxDensity = np.max(lesionVoxels)
                if maxDensity >= 400:
                    densityDistribution[3] += 1
                elif maxDensity >= 300:
                    densityDistribution[2] += 1
                elif maxDensity >= 200:
                    densityDistribution[1] += 1
                else:
                    densityDistribution[0] += 1

        print(f"Slices processed: {totalSlicesProcessed}, Slices filtered (area < {MIN_AREA_MM2} mm²): {totalSlicesFiltered}")

        # Convert density distribution to percentages
        if validLesionCount > 0:
            densityDistribution = [100.0 * x / validLesionCount for x in densityDistribution]

        # Print summary
        print(f"Total lesions detected: {num_lesions}, Valid lesions with area ≥{MIN_AREA_MM2} mm²: {validLesionCount}")

        # If no valid lesions found, return empty results
        if validLesionCount == 0:
            print("No valid calcium lesions detected (all detected regions were smaller than minimum volume threshold)")
            self.currentResults = self.getEmptyResults()
            return self.currentResults

        # Calculate statistics
        meanDensity = np.mean(allDensities) if allDensities else 0
        maxDensity = np.max(allDensities) if allDensities else 0

        # Apply slice normalization factor to score
        # This normalizes the score to 3mm equivalent
        normalizedAgatston = totalAgatston * sliceNormalizationFactor
        print(f"Raw score: {totalAgatston:.1f}, Normalized score: {normalizedAgatston:.1f} (factor: {sliceNormalizationFactor:.3f})")
        totalAgatston = normalizedAgatston

        # Calculate equivalent mass (mg)
        # Standard formula: mass = volume (mm³) × mean_density (HU) × CT_calibration_factor
        # For calcium scoring, typical calibration factor is ~0.001 (1 mg/mm³ at 1000 HU)
        # More accurate: mass = volume × (density_HU / 1000) × tissue_density
        # For hydroxyapatite (calcium): density ≈ 1.2 mg/mm³ at 1000 HU
        equivalentMass_mg = totalVolume * (meanDensity / 1000.0) * 1.2 if meanDensity > 0 else 0

        # Classify severity
        classification, severity = self.classifySeverity(totalAgatston, patientSex)

        # Store results (using validLesionCount instead of num_lesions)
        self.currentResults = {
            'agatston_score': totalAgatston,
            'total_volume_mm3': totalVolume,
            'equivalent_mass_mg': equivalentMass_mg,
            'num_lesions': validLesionCount,
            'mean_density': meanDensity,
            'max_density': maxDensity,
            'density_distribution': densityDistribution,
            'classification': classification,
            'severity': severity,
            'patient_sex': patientSex,
            'patient_age': patientAge,
            'all_densities': allDensities
        }

        return self.currentResults

    def classifySeverity(self, agatstonScore, sex):
        """Classify severity based on literature thresholds (Non-Contrast CT only)

        Args:
            agatstonScore: Calculated Agatston score
            sex: 'M' or 'F' for sex-specific thresholds

        Returns:
            Tuple of (classification_text, severity_level)
        """
        # Non-Contrast CT thresholds based on JAHA 2024 (Tastet et al.)
        # Women ≥1300 AU, Men ≥2000 AU for severe AS

        if sex == 'F':
            if agatstonScore >= 1300:
                return "Severe Aortic Stenosis (AS)", "Severe"
            elif agatstonScore >= 400:
                return "Moderate Aortic Stenosis", "Moderate"
            elif agatstonScore >= 100:
                return "Mild Aortic Valve Calcification", "Mild"
            else:
                return "Minimal/Normal", "Normal/Minimal"
        else:  # Male
            if agatstonScore >= 2000:
                return "Severe Aortic Stenosis (AS)", "Severe"
            elif agatstonScore >= 1000:
                return "Moderate Aortic Stenosis", "Moderate"
            elif agatstonScore >= 100:
                return "Mild Aortic Valve Calcification", "Mild"
            else:
                return "Minimal/Normal", "Normal/Minimal"

    def getEmptyResults(self):
        """Return empty results structure"""
        return {
            'agatston_score': 0,
            'total_volume_mm3': 0,
            'equivalent_mass_mg': 0,
            'num_lesions': 0,
            'mean_density': 0,
            'max_density': 0,
            'density_distribution': [0, 0, 0, 0],
            'classification': "No calcification detected",
            'severity': "Normal/Minimal",
            'patient_sex': 'M',
            'patient_age': None,
            'ct_mode': 'non-contrast',
            'all_densities': []
        }

    def createThoraxSegment(self, segmentationNode, volumeNode):
        """Create a semi-transparent thorax segment for anatomical context in 3D view"""
        try:
            segmentation = segmentationNode.GetSegmentation()

            # Check if thorax segment already exists
            thoraxSegmentId = None
            for i in range(segmentation.GetNumberOfSegments()):
                segmentId = segmentation.GetNthSegmentID(i)
                segment = segmentation.GetSegment(segmentId)
                if segment.GetName() == "Thorax_Context":
                    thoraxSegmentId = segmentId
                    print("✓ Thorax segment already exists, updating...")
                    break

            # Get volume array
            volumeArray = slicer.util.arrayFromVolume(volumeNode)

            # Create thorax mask: soft tissues (between -200 and +200 HU)
            # This captures heart, aorta, mediastinum, and chest wall
            thoraxMask = ((volumeArray >= -200) & (volumeArray <= 200)).astype(np.uint8)

            # Optional: morphological operations to clean up the mask
            from scipy import ndimage
            # Fill small holes
            thoraxMask = ndimage.binary_fill_holes(thoraxMask).astype(np.uint8)
            # Smooth the surface
            thoraxMask = ndimage.binary_closing(thoraxMask, structure=np.ones((3,3,3))).astype(np.uint8)

            # Create or update segment
            if thoraxSegmentId is None:
                thoraxSegmentId = segmentation.AddEmptySegment("Thorax_Context")
                print(f"✓ Created new Thorax segment: {thoraxSegmentId}")

            segment = segmentation.GetSegment(thoraxSegmentId)

            # Set color: light gray/beige for anatomical context
            segment.SetColor(0.85, 0.82, 0.78)  # Beige/bone color

            # Update segment data
            slicer.util.updateSegmentBinaryLabelmapFromArray(
                thoraxMask,
                segmentationNode,
                thoraxSegmentId,
                volumeNode
            )

            # Set transparency for thorax segment
            segmentationNode.CreateClosedSurfaceRepresentation()
            displayNode = segmentationNode.GetDisplayNode()
            if displayNode:
                # Set thorax to 20% opacity (80% transparent)
                displayNode.SetSegmentOpacity3D(thoraxSegmentId, 0.2)
                displayNode.SetVisibility(True)
                displayNode.SetVisibility3D(True)
                print("✓ Thorax segment configured with 20% opacity")

            print(f"✓ Thorax segment created successfully for anatomical context")

        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"Warning: Could not create thorax segment: {str(e)}")
            # Don't fail the entire 3D visualization if thorax creation fails

    def create3DModel(self, segmentationNode, volumeNode, colorByDensity=True):
        """Create 3D model of calcifications with density-based coloring and thorax context"""
        try:
            # Create thorax segment for anatomical context
            self.createThoraxSegment(segmentationNode, volumeNode)

            if colorByDensity and self.currentResults and 'all_densities' in self.currentResults:
                # Split segmentation by density factors
                self.splitSegmentationByDensity(segmentationNode, volumeNode)
            else:
                # Simple single-color 3D display
                segmentationNode.CreateClosedSurfaceRepresentation()
                displayNode = segmentationNode.GetDisplayNode()
                if displayNode:
                    displayNode.SetVisibility(True)
                    displayNode.SetVisibility3D(True)

                    # Set proper opacity for all segments
                    segmentation = segmentationNode.GetSegmentation()
                    for i in range(segmentation.GetNumberOfSegments()):
                        segmentId = segmentation.GetNthSegmentID(i)
                        segment = segmentation.GetSegment(segmentId)
                        segmentName = segment.GetName()

                        # Set opacity based on segment type
                        if "Thorax_Context" in segmentName:
                            displayNode.SetSegmentOpacity3D(segmentId, 0.2)  # 20% for thorax
                        else:
                            displayNode.SetSegmentOpacity3D(segmentId, 1.0)  # 100% for calcium
                            # Set calcium to bright red color
                            segment.SetColor(1.0, 0.0, 0.0)

            # Switch to 3D view layout
            layoutManager = slicer.app.layoutManager()
            layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)

            # Center 3D view on segmentation
            threeDWidget = layoutManager.threeDWidget(0)
            threeDView = threeDWidget.threeDView()
            threeDView.resetFocalPoint()

            slicer.util.infoDisplay("3D model created successfully with thorax context!")

        except Exception as e:
            import traceback
            traceback.print_exc()
            slicer.util.errorDisplay(f"Error creating 3D model:\n{str(e)}")

    def splitSegmentationByDensity(self, segmentationNode, volumeNode):
        """Split calcium segmentation into 4 segments by density factor"""
        try:
            # Get segmentation data
            segmentation = segmentationNode.GetSegmentation()

            # Get the original calcium segment
            originalSegmentId = segmentation.GetNthSegmentID(0)

            # Get volume array and segmentation array
            volumeArray = slicer.util.arrayFromVolume(volumeNode)
            segmentArray = slicer.util.arrayFromSegmentBinaryLabelmap(
                segmentationNode, originalSegmentId, volumeNode
            )

            # Create mask for calcium voxels
            calciumMask = segmentArray > 0

            # Define density ranges and colors
            densityRanges = [
                (130, 199, "Factor1_130-199HU", [1.0, 1.0, 0.0]),    # Yellow
                (200, 299, "Factor2_200-299HU", [1.0, 0.647, 0.0]),  # Orange
                (300, 399, "Factor3_300-399HU", [1.0, 0.4, 0.0]),    # Dark Orange
                (400, 9999, "Factor4_400+HU", [1.0, 0.0, 0.0])       # Red
            ]

            # Hide original segment
            segmentation.GetSegment(originalSegmentId).SetTag("Visibility", "0")

            # Create segments for each density range
            for minHU, maxHU, name, color in densityRanges:
                # Create mask for this density range
                densityMask = calciumMask & (volumeArray >= minHU) & (volumeArray <= maxHU)

                if np.any(densityMask):
                    # Create new segment
                    segmentId = segmentation.AddEmptySegment(name)
                    segment = segmentation.GetSegment(segmentId)
                    segment.SetColor(color[0], color[1], color[2])

                    # Set segment data
                    slicer.util.updateSegmentBinaryLabelmapFromArray(
                        densityMask.astype(np.uint8),
                        segmentationNode,
                        segmentId,
                        volumeNode
                    )

            # Create 3D representation
            segmentationNode.CreateClosedSurfaceRepresentation()

            # Update display
            displayNode = segmentationNode.GetDisplayNode()
            if displayNode:
                displayNode.SetVisibility(True)
                displayNode.SetVisibility3D(True)

                # Hide the original segment in 3D
                displayNode.SetSegmentVisibility3D(originalSegmentId, False)

                # Show all density factor segments with full opacity
                for i in range(segmentation.GetNumberOfSegments()):
                    segmentId = segmentation.GetNthSegmentID(i)
                    segment = segmentation.GetSegment(segmentId)
                    segmentName = segment.GetName()

                    # Set full opacity for calcium segments, keep transparency for thorax
                    if "Thorax_Context" in segmentName:
                        displayNode.SetSegmentOpacity3D(segmentId, 0.2)  # 20% opacity for thorax
                        displayNode.SetSegmentVisibility3D(segmentId, True)
                    elif "Factor" in segmentName:
                        displayNode.SetSegmentOpacity3D(segmentId, 1.0)  # 100% opacity for calcium
                        displayNode.SetSegmentVisibility3D(segmentId, True)

            print(f"Created {segmentation.GetNumberOfSegments() - 1} density-based segments with proper opacity")

        except Exception as e:
            import traceback
            traceback.print_exc()
            slicer.util.errorDisplay(f"Error splitting by density:\n{str(e)}")

    def createCharts(self, volumeNode, segmentationNode):
        """Create comprehensive visualization charts"""
        if not self.currentResults or 'all_densities' not in self.currentResults:
            slicer.util.errorDisplay("No results available. Please calculate Agatston score first.")
            return

        densities = self.currentResults['all_densities']
        if not densities:
            slicer.util.errorDisplay("No calcium detected. Cannot create charts.")
            return

        try:
            # Try to import matplotlib
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                from matplotlib.patches import Wedge, Circle, Rectangle
                import matplotlib.patches as mpatches
                matplotlibAvailable = True
            except ImportError:
                matplotlibAvailable = False
                slicer.util.warningDisplay(
                    "Matplotlib not installed.\n\n"
                    "Only text-based tables will be shown.\n"
                    "Install matplotlib from Settings tab for graphical charts."
                )

            if matplotlibAvailable:
                # Create comprehensive multi-chart visualization
                self.createComprehensiveCharts(segmentationNode, volumeNode)
            else:
                # Fallback to text-based output
                self.createTextBasedSummary()

        except Exception as e:
            import traceback
            traceback.print_exc()
            slicer.util.errorDisplay(f"Error creating charts:\n{str(e)}")

    def createBullsEyeAnalysis(self, volumeNode, segmentationNode):
        """Create bull's eye plot from spatial distribution analysis"""
        if not self.currentResults:
            slicer.util.warningDisplay("Please calculate Agatston score first!")
            return

        try:
            # Try to import matplotlib
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                matplotlibAvailable = True
            except ImportError:
                matplotlibAvailable = False
                slicer.util.warningDisplay(
                    "Matplotlib not installed.\n\n"
                    "Install matplotlib from Settings tab to view bull's eye plot."
                )
                return

            if matplotlibAvailable:
                # Analyze spatial distribution
                regionData = self.analyzeSpatialDistribution(segmentationNode, volumeNode)
                if regionData:
                    # Create bull's eye plot
                    self.createBullsEyePlot(regionData)
                else:
                    slicer.util.warningDisplay("No calcium detected for spatial analysis!")

        except Exception as e:
            import traceback
            traceback.print_exc()
            slicer.util.errorDisplay(f"Error creating bull's eye plot:\n{str(e)}")

    def createComprehensiveCharts(self, segmentationNode, volumeNode):
        """Create a comprehensive 2x2 grid of charts"""
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import Wedge
        import numpy as np
        import tempfile

        # Create 2x2 subplot grid
        fig = plt.figure(figsize=(16, 12))

        # Chart 1: Density Factor Distribution (Top Left)
        ax1 = plt.subplot(2, 2, 1)
        self.createDensityDistributionChart(ax1)

        # Chart 2: Axial Distribution (Top Right)
        ax2 = plt.subplot(2, 2, 2)
        self.createAxialDistributionChart(ax2, segmentationNode, volumeNode)

        # Chart 3: Risk Gauge (Bottom Left)
        ax3 = plt.subplot(2, 2, 3)
        self.createRiskGauge(ax3)

        # Chart 4: Age/Sex Percentile (Bottom Right)
        ax4 = plt.subplot(2, 2, 4)
        self.createPercentileChart(ax4)

        # Overall title
        agatstonScore = self.currentResults.get('agatston_score', 0)
        severity = self.currentResults.get('severity', 'Unknown')
        fig.suptitle(f'Aortic Valve Calcium Score Analysis\nAgatston Score: {agatstonScore:.1f} AU - {severity}',
                    fontsize=16, fontweight='bold')

        plt.tight_layout(rect=[0, 0.03, 1, 0.96])

        # Save to temporary file
        tempDir = tempfile.gettempdir()
        plotPath = os.path.join(tempDir, "calcium_analysis_comprehensive.png")
        plt.savefig(plotPath, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"✓ Comprehensive charts saved to: {plotPath}")

        # Display in Slicer
        plotVolume = slicer.util.loadVolume(plotPath)
        if plotVolume:
            layoutManager = slicer.app.layoutManager()
            layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)
            slicer.util.setSliceViewerLayers(background=plotVolume, fit=True)
            print("✓ Charts displayed in Slicer")
            slicer.util.infoDisplay(
                "Comprehensive charts created successfully!\n\n"
                "The visualization includes:\n"
                "1. Density Factor Distribution\n"
                "2. Axial Distribution (slice-by-slice)\n"
                "3. Risk Stratification Gauge\n"
                "4. Age/Sex Percentile Comparison"
            )

    def createDensityDistributionChart(self, ax):
        """Chart 1: Bar chart showing distribution of lesions by density factor"""
        import numpy as np

        # Get density distribution from results
        densityDist = self.currentResults.get('density_distribution', [0, 0, 0, 0])

        factors = ['Factor 1\n(130-199 HU)', 'Factor 2\n(200-299 HU)',
                   'Factor 3\n(300-399 HU)', 'Factor 4\n(400+ HU)']
        colors = ['#fff44f', '#ffaa00', '#ff6600', '#ff0000']

        bars = ax.bar(factors, densityDist, color=colors, edgecolor='black', linewidth=1.5)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.1f}%',
                       ha='center', va='bottom', fontweight='bold', fontsize=10)

        ax.set_ylabel('Percentage of Lesions (%)', fontsize=12, fontweight='bold')
        ax.set_title('Density Factor Distribution', fontsize=14, fontweight='bold')
        ax.set_ylim(0, max(densityDist) * 1.2 if max(densityDist) > 0 else 100)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

    def createAxialDistributionChart(self, ax, segmentationNode, volumeNode):
        """Chart 2: Line plot showing calcium score distribution along z-axis"""
        import numpy as np

        # Get segment array
        segmentation = segmentationNode.GetSegmentation()
        segmentId = segmentation.GetNthSegmentID(0)
        volumeArray = slicer.util.arrayFromVolume(volumeNode)
        segmentArray = slicer.util.arrayFromSegmentBinaryLabelmap(
            segmentationNode, segmentId, volumeNode
        )

        # Calculate score per slice
        spacing = volumeNode.GetSpacing()
        sliceArea = spacing[0] * spacing[1]

        scorePerSlice = []
        sliceIndices = []

        for sliceIdx in range(segmentArray.shape[0]):
            sliceMask = segmentArray[sliceIdx, :, :]
            if np.any(sliceMask):
                pixelCount = np.sum(sliceMask)
                area_mm2 = pixelCount * sliceArea

                # Get max density in slice
                sliceVoxels = volumeArray[sliceIdx, :, :][sliceMask]
                maxDensity = np.max(sliceVoxels)

                # Density factor
                if maxDensity >= 400:
                    densityFactor = 4
                elif maxDensity >= 300:
                    densityFactor = 3
                elif maxDensity >= 200:
                    densityFactor = 2
                else:
                    densityFactor = 1

                score = area_mm2 * densityFactor
                scorePerSlice.append(score)
                sliceIndices.append(sliceIdx)

        if scorePerSlice:
            ax.plot(sliceIndices, scorePerSlice, marker='o', linewidth=2,
                   markersize=4, color='#0078d7')
            ax.fill_between(sliceIndices, scorePerSlice, alpha=0.3, color='#0078d7')

            ax.set_xlabel('Slice Index (Inferior → Superior)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Agatston Score', fontsize=12, fontweight='bold')
            ax.set_title('Axial Distribution of Calcium', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
        else:
            ax.text(0.5, 0.5, 'No calcium detected', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)

    def createRiskGauge(self, ax):
        """Chart 3: Risk stratification gauge (speedometer style)"""
        import numpy as np
        from matplotlib.patches import Wedge, Circle

        agatstonScore = self.currentResults.get('agatston_score', 0)
        sex = self.currentResults.get('patient_sex', 'Unknown')

        # Severity thresholds based on JAHA 2024
        if sex == 'F' or sex == 'Female':
            thresholds = [0, 400, 1300, 5000]
            labels = ['Minimal/Mild', 'Moderate', 'Severe']
        else:  # Male or Unknown
            thresholds = [0, 1000, 2000, 5000]
            labels = ['Minimal/Mild', 'Moderate', 'Severe']

        colors = ['#28a745', '#ffc107', '#dc3545']

        # Draw gauge arcs
        startAngle = 180
        totalAngle = 180
        anglePerZone = totalAngle / 3

        for i in range(3):
            wedge = Wedge((0, 0), 1.0,
                         startAngle + i * anglePerZone,
                         startAngle + (i + 1) * anglePerZone,
                         width=0.3, facecolor=colors[i],
                         edgecolor='white', linewidth=2)
            ax.add_patch(wedge)

        # Calculate needle angle
        if agatstonScore < thresholds[1]:
            # Minimal/Mild zone
            proportion = agatstonScore / thresholds[1]
            needleAngle = 180 + proportion * anglePerZone
        elif agatstonScore < thresholds[2]:
            # Moderate zone
            proportion = (agatstonScore - thresholds[1]) / (thresholds[2] - thresholds[1])
            needleAngle = 180 + anglePerZone + proportion * anglePerZone
        else:
            # Severe zone
            proportion = min((agatstonScore - thresholds[2]) / (thresholds[3] - thresholds[2]), 1.0)
            needleAngle = 180 + 2 * anglePerZone + proportion * anglePerZone

        # Draw needle
        needleRad = np.radians(needleAngle)
        ax.arrow(0, 0, 0.7 * np.cos(needleRad), 0.7 * np.sin(needleRad),
                head_width=0.1, head_length=0.1, fc='black', ec='black', linewidth=2)

        # Center circle
        center = Circle((0, 0), 0.1, color='black')
        ax.add_patch(center)

        # Labels
        ax.text(0, -0.3, f'{agatstonScore:.0f} AU',
               ha='center', va='top', fontsize=20, fontweight='bold')
        ax.text(0, -0.5, self.currentResults.get('severity', ''),
               ha='center', va='top', fontsize=14, color='#666')

        # Zone labels
        ax.text(-0.85, 0.15, 'Minimal/\nMild', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(0, 0.85, 'Moderate', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(0.85, 0.15, 'Severe', ha='center', va='center', fontsize=9, fontweight='bold')

        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-0.7, 1.3)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title('Risk Stratification', fontsize=14, fontweight='bold', pad=20)

    def createPercentileChart(self, ax):
        """Chart 4: Age/Sex percentile comparison based on MESA study"""
        import numpy as np

        agatstonScore = self.currentResults.get('agatston_score', 0)

        # Get age and sex from results
        age = self.currentResults.get('patient_age', None)
        sex = self.currentResults.get('patient_sex', 'Unknown')

        if age is None or age < 40 or age > 85:
            ax.text(0.5, 0.5, 'Enter patient age (40-85 years)\nin Patient Information\nto show percentile comparison',
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title('Age/Sex Percentile Comparison\n(MESA Study Reference)',
                        fontsize=14, fontweight='bold')
            ax.axis('off')
            return

        # MESA study percentile data (simplified, based on published literature)
        # Percentiles: 25th, 50th, 75th, 90th for age groups
        # This is a simplified version - real data would be more granular

        if sex == 'F' or sex == 'Female':
            # Female percentiles by age group (40-44, 45-54, 55-64, 65-74, 75-84)
            percentile_data = {
                (40, 44): [0, 0, 5, 38],
                (45, 54): [0, 0, 26, 95],
                (55, 64): [0, 4, 78, 250],
                (65, 74): [0, 25, 185, 515],
                (75, 84): [0, 60, 330, 800],
            }
            color = '#e91e63'
            sexLabel = 'Female'
        else:
            # Male percentiles
            percentile_data = {
                (40, 44): [0, 0, 11, 84],
                (45, 54): [0, 3, 54, 203],
                (55, 64): [0, 23, 161, 540],
                (65, 74): [0, 81, 384, 950],
                (75, 84): [0, 155, 610, 1400],
            }
            color = '#2196f3'
            sexLabel = 'Male'

        # Find appropriate age group
        for (minAge, maxAge), percentiles in percentile_data.items():
            if minAge <= age <= maxAge:
                p25, p50, p75, p90 = percentiles
                break
        else:
            # Default to oldest group if age > 84
            p25, p50, p75, p90 = list(percentile_data.values())[-1]

        # Plot percentile bars
        categories = ['25th', '50th', '75th', '90th', 'Patient']
        values = [p25, p50, p75, p90, agatstonScore]
        colors_bars = ['#90caf9', '#64b5f6', '#42a5f5', '#1e88e5', color]

        bars = ax.barh(categories, values, color=colors_bars, edgecolor='black', linewidth=1.5)

        # Add value labels
        for i, (bar, val) in enumerate(zip(bars, values)):
            ax.text(val, bar.get_y() + bar.get_height()/2,
                   f'  {val:.0f}',
                   va='center', ha='left', fontweight='bold', fontsize=10)

        # Determine percentile category
        if agatstonScore <= p25:
            percentile_text = '<25th percentile'
        elif agatstonScore <= p50:
            percentile_text = '25th-50th percentile'
        elif agatstonScore <= p75:
            percentile_text = '50th-75th percentile'
        elif agatstonScore <= p90:
            percentile_text = '75th-90th percentile'
        else:
            percentile_text = '>90th percentile'

        ax.set_xlabel('Agatston Score (AU)', fontsize=12, fontweight='bold')
        ax.set_title(f'Percentile Comparison\n{sexLabel}, Age {age} (MESA Study)\n{percentile_text}',
                    fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

    def createTextBasedSummary(self):
        """Fallback text-based summary when matplotlib not available"""
        message = "=== CALCIUM SCORE ANALYSIS ===\n\n"

        message += f"Agatston Score: {self.currentResults.get('agatston_score', 0):.1f} AU\n"
        message += f"Severity: {self.currentResults.get('severity', 'Unknown')}\n"
        message += f"Total Volume: {self.currentResults.get('total_volume', 0):.1f} mm³\n\n"

        message += "Density Distribution:\n"
        densityDist = self.currentResults.get('density_distribution', [0, 0, 0, 0])
        message += f"  Factor 1 (130-199 HU): {densityDist[0]:.1f}%\n"
        message += f"  Factor 2 (200-299 HU): {densityDist[1]:.1f}%\n"
        message += f"  Factor 3 (300-399 HU): {densityDist[2]:.1f}%\n"
        message += f"  Factor 4 (400+ HU):    {densityDist[3]:.1f}%\n"

        slicer.util.infoDisplay(message)
        print(message)

    def analyzeSpatialDistribution(self, segmentationNode, volumeNode):
        """Analyze spatial distribution of calcium across aortic valve regions"""
        import numpy as np

        # Get segmentation data
        segmentation = segmentationNode.GetSegmentation()
        segmentId = segmentation.GetNthSegmentID(0)

        # Get arrays
        volumeArray = slicer.util.arrayFromVolume(volumeNode)
        segmentArray = slicer.util.arrayFromSegmentBinaryLabelmap(
            segmentationNode, segmentId, volumeNode
        )

        # Find calcium voxels
        calciumMask = segmentArray > 0
        calciumIndices = np.argwhere(calciumMask)

        if len(calciumIndices) == 0:
            return None

        # Calculate centroid (center of calcium mass)
        centroid = np.mean(calciumIndices, axis=0)

        # Define regions based on cylindrical coordinates
        # We'll divide the valve into:
        # - Center: annulus (r < 20% of max radius)
        # - Inner ring: 3 cusps at 120° intervals
        # - Outer ring: 6 regions at 60° intervals

        regionData = {
            'center': {'density': [], 'volume': 0},
            'NC': {'density': [], 'volume': 0},      # Non-coronary cusp
            'RC': {'density': [], 'volume': 0},      # Right coronary cusp
            'LC': {'density': [], 'volume': 0},      # Left coronary cusp
            'NC_basal': {'density': [], 'volume': 0},
            'NC_comm': {'density': [], 'volume': 0},
            'RC_basal': {'density': [], 'volume': 0},
            'RC_comm': {'density': [], 'volume': 0},
            'LC_basal': {'density': [], 'volume': 0},
            'LC_comm': {'density': [], 'volume': 0},
        }

        # Calculate maximum radius
        distances = np.sqrt(np.sum((calciumIndices - centroid) ** 2, axis=1))
        maxRadius = np.max(distances)
        innerThreshold = maxRadius * 0.2
        midThreshold = maxRadius * 0.6

        # Classify each calcium voxel
        for idx in calciumIndices:
            # Get relative position
            relPos = idx - centroid
            distance = np.sqrt(np.sum(relPos ** 2))

            # Calculate angle in axial plane (ignoring superior-inferior)
            angle = np.arctan2(relPos[2], relPos[1])  # j, i components
            angle_deg = np.degrees(angle) % 360

            # Get density value
            density = volumeArray[tuple(idx)]

            # Assign to region based on distance and angle
            if distance < innerThreshold:
                # Center region (annulus)
                regionData['center']['density'].append(density)
                regionData['center']['volume'] += 1
            elif distance < midThreshold:
                # Inner ring - 3 cusps
                if 30 <= angle_deg < 150:
                    regionData['NC']['density'].append(density)
                    regionData['NC']['volume'] += 1
                elif 150 <= angle_deg < 270:
                    regionData['RC']['density'].append(density)
                    regionData['RC']['volume'] += 1
                else:  # 270-30 (wrapping around)
                    regionData['LC']['density'].append(density)
                    regionData['LC']['volume'] += 1
            else:
                # Outer ring - 6 regions
                if 0 <= angle_deg < 60:
                    regionData['LC_comm']['density'].append(density)
                    regionData['LC_comm']['volume'] += 1
                elif 60 <= angle_deg < 120:
                    regionData['NC_basal']['density'].append(density)
                    regionData['NC_basal']['volume'] += 1
                elif 120 <= angle_deg < 180:
                    regionData['NC_comm']['density'].append(density)
                    regionData['NC_comm']['volume'] += 1
                elif 180 <= angle_deg < 240:
                    regionData['RC_basal']['density'].append(density)
                    regionData['RC_basal']['volume'] += 1
                elif 240 <= angle_deg < 300:
                    regionData['RC_comm']['density'].append(density)
                    regionData['RC_comm']['volume'] += 1
                else:  # 300-360
                    regionData['LC_basal']['density'].append(density)
                    regionData['LC_basal']['volume'] += 1

        # Calculate mean density for each region
        for region in regionData.keys():
            if regionData[region]['density']:
                regionData[region]['mean_density'] = np.mean(regionData[region]['density'])
            else:
                regionData[region]['mean_density'] = 0

        print("✓ Spatial distribution analysis completed")
        return regionData

    def createBullsEyePlot(self, regionData):
        """Create bull's eye plot using matplotlib"""
        # Try to import matplotlib, if not available provide instructions
        try:
            import matplotlib
            matplotlib.use('Agg')  # Use non-interactive backend
            import matplotlib.pyplot as plt
            from matplotlib.patches import Wedge, Circle
        except ImportError:
            slicer.util.errorDisplay(
                "Matplotlib is not installed in Slicer.\n\n"
                "To install matplotlib:\n"
                "1. Open the Python Interactor (View > Python Interactor)\n"
                "2. Run this command:\n"
                "   pip_install('matplotlib')\n\n"
                "Or create a text-based distribution table instead?"
            )
            # Fallback: create simple text table
            self.createDistributionTable(regionData)
            return

        import tempfile

        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw=dict(aspect="equal"))

        # Color map based on mean density
        def get_color(mean_density):
            if mean_density == 0:
                return '#f0f0f0'  # Gray for no calcium
            elif mean_density < 200:
                return '#fff44f'  # Yellow - Factor 1
            elif mean_density < 300:
                return '#ffaa00'  # Orange - Factor 2
            elif mean_density < 400:
                return '#ff6600'  # Dark Orange - Factor 3
            else:
                return '#ff0000'  # Red - Factor 4

        # Draw center circle (annulus)
        center_color = get_color(regionData['center']['mean_density'])
        center_circle = Circle((0, 0), 0.2, facecolor=center_color, edgecolor='white', linewidth=2)
        ax.add_patch(center_circle)

        # Inner ring - 3 cusps (120° each)
        cusps = [
            ('NC', 30, 150),   # Non-coronary: 30-150°
            ('RC', 150, 270),  # Right coronary: 150-270°
            ('LC', 270, 390),  # Left coronary: 270-390° (wraps to 30°)
        ]

        for cusp_name, start_angle, end_angle in cusps:
            color = get_color(regionData[cusp_name]['mean_density'])
            wedge = Wedge((0, 0), 0.6, start_angle, end_angle, width=0.4,
                         facecolor=color, edgecolor='white', linewidth=2)
            ax.add_patch(wedge)

        # Outer ring - 6 regions (60° each)
        outer_regions = [
            ('LC_comm', 0, 60),
            ('NC_basal', 60, 120),
            ('NC_comm', 120, 180),
            ('RC_basal', 180, 240),
            ('RC_comm', 240, 300),
            ('LC_basal', 300, 360),
        ]

        for region_name, start_angle, end_angle in outer_regions:
            color = get_color(regionData[region_name]['mean_density'])
            wedge = Wedge((0, 0), 1.0, start_angle, end_angle, width=0.4,
                         facecolor=color, edgecolor='white', linewidth=2)
            ax.add_patch(wedge)

        # Add labels
        ax.text(0, 0, 'Annulus', ha='center', va='center', fontsize=10, weight='bold')

        # Cusp labels (inner ring)
        ax.text(0, 0.4, 'NC', ha='center', va='center', fontsize=12, weight='bold')
        ax.text(-0.35, -0.2, 'RC', ha='center', va='center', fontsize=12, weight='bold')
        ax.text(0.35, -0.2, 'LC', ha='center', va='center', fontsize=12, weight='bold')

        # Add title and legend
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.axis('off')
        ax.set_title('Aortic Valve Calcium Distribution - Bull\'s Eye View',
                    fontsize=16, weight='bold', pad=20)

        # Add colorbar legend
        from matplotlib.patches import Rectangle
        legend_elements = [
            Rectangle((0, 0), 1, 1, fc='#f0f0f0', label='No calcium'),
            Rectangle((0, 0), 1, 1, fc='#fff44f', label='Factor 1 (130-199 HU)'),
            Rectangle((0, 0), 1, 1, fc='#ffaa00', label='Factor 2 (200-299 HU)'),
            Rectangle((0, 0), 1, 1, fc='#ff6600', label='Factor 3 (300-399 HU)'),
            Rectangle((0, 0), 1, 1, fc='#ff0000', label='Factor 4 (400+ HU)'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=10)

        # Save to temporary file
        tempDir = tempfile.gettempdir()
        plotPath = os.path.join(tempDir, "aortic_valve_bullseye.png")
        plt.tight_layout()
        plt.savefig(plotPath, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"✓ Bull's eye plot saved to: {plotPath}")

        # Display in Slicer
        # Load as scalar volume to display in slice views
        plotVolume = slicer.util.loadVolume(plotPath)
        if plotVolume:
            # Switch to layout that shows the image
            layoutManager = slicer.app.layoutManager()
            layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)

            # Set the plot as the background
            slicer.util.setSliceViewerLayers(background=plotVolume, fit=True)

            print("✓ Bull's eye plot displayed in Slicer")
        else:
            # Fallback: just show a message with the path
            slicer.util.infoDisplay(f"Bull's eye plot created at:\n{plotPath}\n\nYou can open it with an external viewer.")

    def createDistributionTable(self, regionData):
        """Create a text-based table showing calcium distribution (fallback when matplotlib unavailable)"""
        # Create a simple table showing distribution
        message = "Aortic Valve Calcium Distribution:\n\n"
        message += "=" * 60 + "\n"
        message += f"{'Region':<25} {'Mean HU':<12} {'Volume (voxels)':<15}\n"
        message += "=" * 60 + "\n"

        # Group regions
        message += "\nCENTER:\n"
        message += f"  Annulus               {regionData['center']['mean_density']:>8.1f} HU    {regionData['center']['volume']:>8} voxels\n"

        message += "\nINNER RING (Cusps):\n"
        message += f"  Non-Coronary (NC)     {regionData['NC']['mean_density']:>8.1f} HU    {regionData['NC']['volume']:>8} voxels\n"
        message += f"  Right Coronary (RC)   {regionData['RC']['mean_density']:>8.1f} HU    {regionData['RC']['volume']:>8} voxels\n"
        message += f"  Left Coronary (LC)    {regionData['LC']['mean_density']:>8.1f} HU    {regionData['LC']['volume']:>8} voxels\n"

        message += "\nOUTER RING (Regions):\n"
        message += f"  NC Basal              {regionData['NC_basal']['mean_density']:>8.1f} HU    {regionData['NC_basal']['volume']:>8} voxels\n"
        message += f"  NC Commissural        {regionData['NC_comm']['mean_density']:>8.1f} HU    {regionData['NC_comm']['volume']:>8} voxels\n"
        message += f"  RC Basal              {regionData['RC_basal']['mean_density']:>8.1f} HU    {regionData['RC_basal']['volume']:>8} voxels\n"
        message += f"  RC Commissural        {regionData['RC_comm']['mean_density']:>8.1f} HU    {regionData['RC_comm']['volume']:>8} voxels\n"
        message += f"  LC Basal              {regionData['LC_basal']['mean_density']:>8.1f} HU    {regionData['LC_basal']['volume']:>8} voxels\n"
        message += f"  LC Commissural        {regionData['LC_comm']['mean_density']:>8.1f} HU    {regionData['LC_comm']['volume']:>8} voxels\n"

        message += "\n" + "=" * 60 + "\n"
        message += "\nDensity Legend:\n"
        message += "  < 200 HU  = Factor 1 (Yellow)\n"
        message += "  200-299   = Factor 2 (Orange)\n"
        message += "  300-399   = Factor 3 (Dark Orange)\n"
        message += "  400+      = Factor 4 (Red)\n"

        # Show in a dialog
        slicer.util.infoDisplay(message)
        print(message)

    def generateAxialMIP(self, volumeNode, segmentationNode, outputPath, timestamp):
        """Generate Maximum Intensity Projection (MIP) of axial view with calcium overlay

        Args:
            volumeNode: CT volume node
            segmentationNode: Segmentation node containing calcium
            outputPath: Directory to save the MIP image
            timestamp: Timestamp string for filename

        Returns:
            str: Path to saved MIP image, or None if failed
        """
        try:
            import numpy as np
            import tempfile

            # Get volume array
            volumeArray = slicer.util.arrayFromVolume(volumeNode)

            # Get segmentation array
            segmentId = segmentationNode.GetSegmentation().GetNthSegmentID(0)
            segmentArray = slicer.util.arrayFromSegmentBinaryLabelmap(
                segmentationNode, segmentId, volumeNode
            )

            # Find slices with calcium
            slicesWithCalcium = np.any(segmentArray, axis=(1, 2))
            calciumSliceIndices = np.where(slicesWithCalcium)[0]

            if len(calciumSliceIndices) == 0:
                print("No calcium found for MIP generation")
                return None

            # Calculate MIP range (expand by 10 slices on each side for context)
            minSlice = max(0, calciumSliceIndices[0] - 10)
            maxSlice = min(volumeArray.shape[0], calciumSliceIndices[-1] + 11)

            print(f"Generating MIP from slices {minSlice} to {maxSlice} (calcium range: {calciumSliceIndices[0]}-{calciumSliceIndices[-1]})")

            # Create MIP (Maximum Intensity Projection)
            mipImage = np.max(volumeArray[minSlice:maxSlice, :, :], axis=0)

            # Create calcium overlay (also MIP of segmentation)
            calciumMIP = np.max(segmentArray[minSlice:maxSlice, :, :], axis=0)

            # Apply window/level for visualization (calcium preset)
            windowLevel = 300
            windowWidth = 1500
            minHU = windowLevel - windowWidth / 2
            maxHU = windowLevel + windowWidth / 2

            # Normalize to 0-255 range
            mipNormalized = np.clip((mipImage - minHU) / (maxHU - minHU) * 255, 0, 255).astype(np.uint8)

            # Create RGB image
            from PIL import Image as PILImage

            # Convert to RGB
            mipRGB = np.stack([mipNormalized, mipNormalized, mipNormalized], axis=-1)

            # Overlay calcium in red
            mipRGB[calciumMIP > 0, 0] = 255  # Red channel
            mipRGB[calciumMIP > 0, 1] = np.clip(mipRGB[calciumMIP > 0, 1] * 0.3, 0, 255)  # Reduce green
            mipRGB[calciumMIP > 0, 2] = np.clip(mipRGB[calciumMIP > 0, 2] * 0.3, 0, 255)  # Reduce blue

            # Create PIL Image (note: need to flip for correct orientation)
            pilImage = PILImage.fromarray(mipRGB)

            # Rotate 180 degrees for correct orientation in PDF
            pilImage = pilImage.rotate(180)

            # Save to file
            mipPath = os.path.join(outputPath, f"AxialMIP_{timestamp}.png")
            pilImage.save(mipPath)

            print(f"✓ Axial MIP saved to: {mipPath}")
            return mipPath

        except Exception as e:
            print(f"Error generating axial MIP: {str(e)}")
            import traceback
            traceback.print_exc()
            return None

    def generatePDFReport(self, volumeNode, segmentationNode, patientInfo, outputPath, options):
        """Generate PDF report"""
        try:
            # Check if results are available
            if not self.currentResults or 'severity' not in self.currentResults:
                slicer.util.errorDisplay(
                    "No results available to generate report.\n\n"
                    "Please calculate the Agatston Score first (Tab 3)."
                )
                return False

            from reportlab.lib import colors
            from reportlab.lib.pagesizes import letter, A4
            from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image, PageBreak
            from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
            from reportlab.lib.units import inch
            from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY

            reportPath = os.path.join(outputPath, f"AorticValveCalcium_Report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf")

            # Create PDF document
            doc = SimpleDocTemplate(reportPath, pagesize=letter,
                                  rightMargin=72, leftMargin=72,
                                  topMargin=72, bottomMargin=18)

            # Container for PDF elements
            story = []
            styles = getSampleStyleSheet()

            # Custom styles
            titleStyle = ParagraphStyle(
                'CustomTitle',
                parent=styles['Heading1'],
                fontSize=20,
                textColor=colors.black,
                spaceAfter=20,
                alignment=TA_CENTER,
                fontName='Helvetica-Bold'
            )

            headingStyle = ParagraphStyle(
                'CustomHeading',
                parent=styles['Heading2'],
                fontSize=12,
                textColor=colors.HexColor('#0078d7'),
                spaceAfter=8,
                spaceBefore=8,
                fontName='Helvetica-Bold'
            )

            # Logo and company info (placeholder - will add to settings later)
            logoStyle = ParagraphStyle(
                'Logo',
                parent=styles['Normal'],
                fontSize=9,
                textColor=colors.grey,
                alignment=TA_LEFT,
                fontName='Helvetica-Oblique'
            )
            story.append(Paragraph("LOGO AZIENDA DA INSERIRE NEI SETTING", logoStyle))
            story.append(Paragraph("Descrizione Azienda", logoStyle))
            story.append(Spacer(1, 0.15*inch))

            # Title
            story.append(Paragraph("AORTIC VALVE – CALCIUM SCORE", titleStyle))
            story.append(Spacer(1, 0.15*inch))

            # Patient Information Section
            story.append(Paragraph("PATIENT INFORMATION", headingStyle))

            # Analysis ID
            analysisID = patientInfo.get('id', 'N/A')
            story.append(Paragraph(f"Analysis ID <b>[{analysisID}]</b>", styles['Normal']))
            story.append(Spacer(1, 0.05*inch))

            # Patient info table - single row with Name, Sex, Age
            patientData = [
                ['Name', 'Sex', 'Age'],
                [patientInfo['name'], patientInfo['sex'], patientInfo['age']]
            ]

            patientTable = Table(patientData, colWidths=[3*inch, 1.5*inch, 1*inch])
            patientTable.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.white),
                ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTNAME', (0, 1), (-1, 1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                ('TOPPADDING', (0, 0), (-1, -1), 4),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
                ('LINEBELOW', (0, 0), (-1, 0), 1, colors.black)
            ]))
            story.append(patientTable)
            story.append(Spacer(1, 0.08*inch))

            # Analysis date
            story.append(Paragraph(f"Analysis ID    <b>[{analysisID}]</b>", styles['Normal']))
            story.append(Paragraph(f"Analysis date  <b>[{patientInfo['date']}]</b>", styles['Normal']))
            story.append(Spacer(1, 0.15*inch))

            # Classification Section
            results = self.currentResults
            story.append(Paragraph("CLASSIFICATION", headingStyle))

            # Set severity color based on result
            severityColor = colors.green
            if results['severity'] == 'Severe':
                severityColor = colors.red
            elif results['severity'] == 'Moderate':
                severityColor = colors.orange
            elif results['severity'] == 'Mild':
                severityColor = colors.yellow

            classificationData = [
                ['Category', 'Assessment'],
                ['Severity', f"{results['severity'].upper()} *"],
                ['Classification', f"{results['classification']} *"]
            ]

            classificationTable = Table(classificationData, colWidths=[2*inch, 4*inch])
            classificationTable.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0078d7')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('BACKGROUND', (1, 1), (1, 1), severityColor),
                ('TEXTCOLOR', (1, 1), (1, 1), colors.white if results['severity'] == 'Severe' else colors.black),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, -1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
                ('TOPPADDING', (0, 0), (-1, -1), 5),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.black)
            ]))
            story.append(classificationTable)
            story.append(Spacer(1, 0.15*inch))

            # Quantitative Results Section
            story.append(Paragraph("QUANTITATIVE RESULTS", headingStyle))

            # CT mode and scoring method (compact format)
            story.append(Paragraph(f"CT Acquisition Mode      <i>[ECG Gated Non-Contrast CT]</i>", styles['Normal']))
            story.append(Paragraph(f"Scoring Method           <i>Standard Agatston (130 HU)</i>", styles['Normal']))
            story.append(Spacer(1, 0.08*inch))

            # Main results table - compact 2-column format
            resultsData = [
                ['Agatston Score', f"{results['agatston_score']:.1f} AU"],
                ['Total Calcium volume (mm³)', f"{results['total_volume_mm3']:.2f} mm³"],
                ['Mean Density (HU)', f"{results['mean_density']:.1f} HU"],
                ['Max Density (HU)', f"{results['max_density']:.1f} HU"],
                ['Equivalent Mass (mg)', f"{results['equivalent_mass_mg']:.2f} mg"]
            ]

            resultsTable = Table(resultsData, colWidths=[3*inch, 2.5*inch])
            resultsTable.setStyle(TableStyle([
                ('ALIGN', (0, 0), (0, -1), 'LEFT'),
                ('ALIGN', (1, 0), (1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                ('FONTNAME', (1, 0), (1, -1), 'Helvetica'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                ('TOPPADDING', (0, 0), (-1, -1), 4),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
                ('LINEABOVE', (0, 0), (-1, 0), 1, colors.black)
            ]))
            story.append(resultsTable)
            story.append(Spacer(1, 0.15*inch))

            # Density Distribution Section - compact list format
            story.append(Paragraph("DENSITY DISTRIBUTION", headingStyle))

            densityText = f"Factor 1  130-199 HU  -  <b>{results['density_distribution'][0]:.1f} %</b><br/>"
            densityText += f"Factor 2  200-299 HU  -  <b>{results['density_distribution'][1]:.1f} %</b><br/>"
            densityText += f"Factor 3  300-399 HU  -  <b>{results['density_distribution'][2]:.1f} %</b><br/>"
            densityText += f"Factor 4  400 HU       -  <b>{results['density_distribution'][3]:.1f} %</b>"

            story.append(Paragraph(densityText, styles['Normal']))
            story.append(Spacer(1, 0.15*inch))

            # References Section - compact format with asterisk notation
            story.append(Paragraph("*REFERENCES", headingStyle))
            referencesStyle = ParagraphStyle(
                'References',
                parent=styles['BodyText'],
                fontSize=8,
                spaceAfter=3,
                leftIndent=10
            )
            story.append(Paragraph(" 1. Clavel MA, et al. Impact of Aortic Valve Calcification. Circulation. 2015", referencesStyle))
            story.append(Paragraph(" 2. ACC/AHA Guidelines for the Management of Valvular Heart Disease. 2020.", referencesStyle))
            story.append(Paragraph("* 3. Sex-specific thresholds: Women >1300 AU, Men >2000 AU for severe AS (JAHA 2024).", referencesStyle))
            story.append(Paragraph(" 4. Threshold: 130 HU for non-contrast CT calcium detection.", referencesStyle))
            story.append(Spacer(1, 0.15*inch))

            # Footer - page 1
            footerStyle = ParagraphStyle(
                'Footer',
                parent=styles['Normal'],
                fontSize=8,
                textColor=colors.grey,
                alignment=TA_CENTER
            )
            story.append(Spacer(1, 0.1*inch))
            story.append(Paragraph(f"REPORT GENERATED ON [{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]", footerStyle))
            story.append(Paragraph("Generated with 3D Slicer – <i>Aortic Valve Calcium Score Plugin</i> (Developed V. Censullo – 2025)", footerStyle))

            # PAGE BREAK - All images on subsequent pages
            story.append(PageBreak())
            story.append(Paragraph("ATTACHMENTS – IMAGING VIEWS", headingStyle))

            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

            # Capture 3D view
            try:
                view3DPath = os.path.join(outputPath, f"3DView_{timestamp}.png")
                print(f"Attempting to capture 3D view to: {view3DPath}")

                if self.captureView(view3DPath, '3D'):
                    if os.path.exists(view3DPath):
                        print(f"3D view file exists, size: {os.path.getsize(view3DPath)} bytes")
                        img = Image(view3DPath, width=4.5*inch, height=3.4*inch)
                        story.append(Paragraph("<b>3D Visualization of Calcium Segmentation</b>", styles['Normal']))
                        story.append(Spacer(1, 0.1*inch))
                        story.append(img)
                        story.append(Spacer(1, 0.2*inch))
                    else:
                        print(f"3D view file NOT found at {view3DPath}")
                        story.append(Paragraph("(3D view screenshot could not be captured)", styles['Normal']))
                else:
                    print("captureView returned False for 3D")
                    story.append(Paragraph("(3D view screenshot could not be captured)", styles['Normal']))
            except Exception as e:
                print(f"Error capturing 3D view: {str(e)}")
                import traceback
                traceback.print_exc()
                story.append(Paragraph(f"(3D view error: {str(e)})", styles['Normal']))

            # Generate Axial MIP (Maximum Intensity Projection)
            try:
                print("Generating Axial MIP...")
                axialMIPPath = self.generateAxialMIP(volumeNode, segmentationNode, outputPath, timestamp)

                if axialMIPPath and os.path.exists(axialMIPPath):
                    print(f"Axial MIP file exists, size: {os.path.getsize(axialMIPPath)} bytes")
                    img = Image(axialMIPPath, width=4.5*inch, height=3.4*inch)
                    story.append(Paragraph("<b>Axial Maximum Intensity Projection (MIP) with Calcium Overlay</b>", styles['Normal']))
                    story.append(Spacer(1, 0.1*inch))
                    story.append(img)
                    story.append(Spacer(1, 0.2*inch))
                else:
                    print(f"Axial MIP could not be generated")
                    story.append(Paragraph("(Axial MIP could not be generated)", styles['Normal']))
            except Exception as e:
                print(f"Error generating Axial MIP: {str(e)}")
                import traceback
                traceback.print_exc()
                story.append(Paragraph(f"(Axial MIP error: {str(e)})", styles['Normal']))

            # Generate and include charts if requested
            if options.get('includeCharts', False):
                try:
                    print("Generating comprehensive charts for PDF...")
                    import tempfile

                    # Check if matplotlib is available
                    try:
                        import matplotlib
                        matplotlib.use('Agg')
                        import matplotlib.pyplot as plt
                        matplotlibAvailable = True
                    except ImportError:
                        matplotlibAvailable = False
                        print("Matplotlib not available, skipping charts")

                    if matplotlibAvailable:
                        # Generate comprehensive charts (4-panel)
                        comprehensiveChartsPath = os.path.join(outputPath, f"ComprehensiveCharts_{timestamp}.png")

                        # Create the charts directly to file
                        fig = plt.figure(figsize=(16, 12))

                        # Chart 1: Density Factor Distribution
                        ax1 = plt.subplot(2, 2, 1)
                        self.createDensityDistributionChart(ax1)

                        # Chart 2: Axial Distribution
                        ax2 = plt.subplot(2, 2, 2)
                        self.createAxialDistributionChart(ax2, segmentationNode, volumeNode)

                        # Chart 3: Risk Gauge
                        ax3 = plt.subplot(2, 2, 3)
                        self.createRiskGauge(ax3)

                        # Chart 4: Percentile
                        ax4 = plt.subplot(2, 2, 4)
                        self.createPercentileChart(ax4)

                        # Overall title
                        agatstonScore = self.currentResults.get('agatston_score', 0)
                        severity = self.currentResults.get('severity', 'Unknown')
                        fig.suptitle(f'Aortic Valve Calcium Score Analysis\nAgatston Score: {agatstonScore:.1f} AU - {severity}',
                                    fontsize=16, fontweight='bold')

                        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
                        plt.savefig(comprehensiveChartsPath, dpi=150, bbox_inches='tight')
                        plt.close()

                        # Add to PDF
                        if os.path.exists(comprehensiveChartsPath):
                            story.append(PageBreak())
                            story.append(Paragraph("COMPREHENSIVE ANALYSIS CHARTS", headingStyle))
                            img = Image(comprehensiveChartsPath, width=6.5*inch, height=4.9*inch)
                            story.append(img)
                            story.append(Spacer(1, 0.3*inch))
                            print(f"✓ Comprehensive charts added to PDF")

                        # Generate Bull's Eye plot
                        bullsEyePath = os.path.join(outputPath, f"BullsEye_{timestamp}.png")

                        # Analyze spatial distribution
                        regionData = self.analyzeSpatialDistribution(segmentationNode, volumeNode)

                        if regionData:
                            # Create bull's eye plot directly to file
                            fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw=dict(aspect="equal"))

                            # Color map
                            def get_color(mean_density):
                                if mean_density == 0:
                                    return '#f0f0f0'
                                elif mean_density < 200:
                                    return '#fff44f'
                                elif mean_density < 300:
                                    return '#ffaa00'
                                elif mean_density < 400:
                                    return '#ff6600'
                                else:
                                    return '#ff0000'

                            from matplotlib.patches import Wedge, Circle, Rectangle

                            # Draw center circle
                            center_color = get_color(regionData['center']['mean_density'])
                            center_circle = Circle((0, 0), 0.2, facecolor=center_color, edgecolor='white', linewidth=2)
                            ax.add_patch(center_circle)

                            # Inner ring - 3 cusps
                            cusps = [
                                ('NC', 30, 150),
                                ('RC', 150, 270),
                                ('LC', 270, 390),
                            ]

                            for cusp_name, start_angle, end_angle in cusps:
                                color = get_color(regionData[cusp_name]['mean_density'])
                                wedge = Wedge((0, 0), 0.6, start_angle, end_angle, width=0.4,
                                             facecolor=color, edgecolor='white', linewidth=2)
                                ax.add_patch(wedge)

                            # Outer ring - 6 regions
                            outer_regions = [
                                ('LC_comm', 0, 60),
                                ('NC_basal', 60, 120),
                                ('NC_comm', 120, 180),
                                ('RC_basal', 180, 240),
                                ('RC_comm', 240, 300),
                                ('LC_basal', 300, 360),
                            ]

                            for region_name, start_angle, end_angle in outer_regions:
                                color = get_color(regionData[region_name]['mean_density'])
                                wedge = Wedge((0, 0), 1.0, start_angle, end_angle, width=0.4,
                                             facecolor=color, edgecolor='white', linewidth=2)
                                ax.add_patch(wedge)

                            # Labels
                            ax.text(0, 0, 'Annulus', ha='center', va='center', fontsize=10, weight='bold')
                            ax.text(0, 0.4, 'NC', ha='center', va='center', fontsize=12, weight='bold')
                            ax.text(-0.35, -0.2, 'RC', ha='center', va='center', fontsize=12, weight='bold')
                            ax.text(0.35, -0.2, 'LC', ha='center', va='center', fontsize=12, weight='bold')

                            # Title and legend
                            ax.set_xlim(-1.3, 1.3)
                            ax.set_ylim(-1.3, 1.3)
                            ax.axis('off')
                            ax.set_title('Aortic Valve Calcium Distribution - Bull\'s Eye View',
                                        fontsize=16, weight='bold', pad=20)

                            legend_elements = [
                                Rectangle((0, 0), 1, 1, fc='#f0f0f0', label='No calcium'),
                                Rectangle((0, 0), 1, 1, fc='#fff44f', label='Factor 1 (130-199 HU)'),
                                Rectangle((0, 0), 1, 1, fc='#ffaa00', label='Factor 2 (200-299 HU)'),
                                Rectangle((0, 0), 1, 1, fc='#ff6600', label='Factor 3 (300-399 HU)'),
                                Rectangle((0, 0), 1, 1, fc='#ff0000', label='Factor 4 (400+ HU)'),
                            ]
                            ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=10)

                            plt.tight_layout()
                            plt.savefig(bullsEyePath, dpi=150, bbox_inches='tight')
                            plt.close()

                            # Add to PDF
                            if os.path.exists(bullsEyePath):
                                story.append(Paragraph("SPATIAL DISTRIBUTION - BULL'S EYE VIEW", headingStyle))
                                img = Image(bullsEyePath, width=5.5*inch, height=5.5*inch)
                                story.append(img)
                                story.append(Spacer(1, 0.3*inch))
                                print(f"✓ Bull's eye plot added to PDF")

                except Exception as e:
                    print(f"Error generating charts for PDF: {str(e)}")
                    import traceback
                    traceback.print_exc()
                    story.append(Paragraph(f"(Charts could not be generated: {str(e)})", styles['Normal']))

            # No additional references or footer at end (already on page 1)

            # Build PDF
            doc.build(story)

            print(f"PDF report saved to: {reportPath}")
            slicer.util.infoDisplay(f"PDF report successfully generated!\n\nSaved to:\n{reportPath}")
            return True

        except ImportError as ie:
            print(f"reportlab not installed: {str(ie)}")
            slicer.util.errorDisplay(
                "reportlab library is required for PDF generation.\n\n"
                "Please install it using:\n"
                "pip install reportlab\n\n"
                "Or use Python Slicer's pip:\n"
                "PythonSlicer -m pip install reportlab"
            )
            return False

        except Exception as e:
            print(f"Error generating PDF report: {str(e)}")
            import traceback
            traceback.print_exc()
            slicer.util.errorDisplay(f"Error generating PDF report:\n{str(e)}")
            return False

    def captureView(self, filename, viewName='3D'):
        """Capture screenshot of specified view"""
        try:
            import qt

            # Get the view widget
            layoutManager = slicer.app.layoutManager()

            if viewName == '3D':
                # Capture 3D view
                threeDWidget = layoutManager.threeDWidget(0)
                threeDView = threeDWidget.threeDView()
                threeDView.forceRender()

                # Use the renderWindow to get a QImage
                renderWindow = threeDView.renderWindow()
                renderWindow.Render()

                # Capture as QPixmap (using modern grab() method)
                screenshot = threeDView.grab()
                screenshot.save(filename)
            else:
                # For slice views
                sliceWidget = layoutManager.sliceWidget(viewName)
                sliceView = sliceWidget.sliceView()
                sliceView.forceRender()

                # Capture slice view as QPixmap (using modern grab() method)
                screenshot = sliceView.grab()
                screenshot.save(filename)

            print(f"Screenshot saved to: {filename}")
            return True

        except Exception as e:
            import traceback
            print(f"Error capturing view {viewName}: {str(e)}")
            traceback.print_exc()
            return False

    def captureChartView(self, filename):
        """Capture the current chart/plot view"""
        try:
            import qt

            layoutManager = slicer.app.layoutManager()
            plotWidget = layoutManager.plotWidget(0)

            if plotWidget:
                plotWidget.repaint()
                slicer.app.processEvents()

                # Capture plot widget as QPixmap (using modern grab() method)
                screenshot = plotWidget.grab()
                screenshot.save(filename)
                print(f"Chart screenshot saved to: {filename}")
                return True
            else:
                print("No plot widget available")
                return False

        except Exception as e:
            import traceback
            print(f"Error capturing chart: {str(e)}")
            traceback.print_exc()
            return False

    def getInterpretation(self, results):
        """Generate interpretation text"""
        score = results['agatston_score']
        severity = results['severity']

        interpretation = f"The aortic valve calcium score of {score:.1f} AU "

        if severity == "Severe":
            interpretation += "indicates severe aortic valve calcification, "
            interpretation += "which is highly suggestive of severe aortic stenosis. "
            interpretation += "Clinical correlation with echocardiography is recommended. "
            interpretation += "Patient may be candidate for aortic valve intervention (TAVR/SAVR)."
        elif severity == "Moderate":
            interpretation += "indicates moderate aortic valve calcification. "
            interpretation += "This suggests at least moderate aortic stenosis. "
            interpretation += "Follow-up with echocardiography is recommended for hemodynamic assessment."
        elif severity == "Mild":
            interpretation += "indicates mild aortic valve calcification. "
            interpretation += "This may represent early aortic valve disease. "
            interpretation += "Regular monitoring and risk factor management are recommended."
        else:
            interpretation += "indicates minimal or no significant aortic valve calcification. "
            interpretation += "This is within normal limits."

        return interpretation

#
# SlicerCaScoreTest
#

class AoCaScoreTest(ScriptedLoadableModuleTest):
    """Test cases for the module"""

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_Basic()

    def test_Basic(self):
        """Basic test"""
        self.delayDisplay("Starting basic test")

        logic = AoCaScoreLogic()
        self.assertIsNotNone(logic)

        self.delayDisplay("Test passed")
