# AoCalciumScore - Aortic Valve Calcium Scoring Plugin

3D Slicer extension for automated quantification of aortic valve calcification using the Agatston scoring method.

![Version](https://img.shields.io/badge/version-3.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![3D Slicer](https://img.shields.io/badge/3D%20Slicer-5.0+-red)

## Authors

**Vittorio Censullo**  
**AITeRTC** - *Associazione Italiana Tecnici di Radiologia e Radiodiagnostica e Tomografia Computerizzata ETS*  
2025

## 🌟 Key Features

### Core Functionality
- **Guided Tab-Based Workflow**: Intuitive step-by-step process for Aortic Valve calcium scoring
- **Advanced Segmentation Methods**:
  - **Point and Click**: One-click calcium chunk/Island segmentation (Region Growing) (recommended)
  - **ROI Box**: Rectangular Box threshold based segmentation
  - **Paint Mode**: Manual painting with threshold-guided brush
- **Intelligent Noise Filtering**: Automatic removal of lesions < 1.0 mm2 to exclude artifacts
- **Sex-Specific Classification**: Based on recent literature (JAHA 2024)
  - Women: Severe ≥1377 AU, Moderate 400-1377 AU
  - Men: Severe ≥2062 AU, Moderate 1000-2062 AU
- **Real-Time 3D Visualization**: Color-coded calcium density (4 levels: 130-199, 200-299, 300-399, ≥400 HU)

### Settings & Persistence
- **Configurable Threshold**: Default 130 HU for non-contrast CT (customisable in settings)
- **Persistent Settings**: Logo, company description, and threshold saved across sessions
- **Easy Dependency Management**: One-click installation of required libraries

## 📦 Installation

### Method 1: Manual Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/vcensullo/AoCalciumScore.git
   ```

2. Open 3D Slicer

3. Go to **Edit → Application Settings → Modules → Additional module paths**

4. Add the path to the `AoCaScore` folder

5. Restart 3D Slicer

6. The module will appear under **Modules → Cardiac → Aortic Valve Calcium Score**


## 🔧 Dependencies

The plugin requires the `reportlab` and 'matplotlib' libraries for PDF and plot generation.

1. Open the plugin in 3D Slicer
2. Go to **Tab 0 - Settings**
3. Click **Install/Update Dependencies**


## 📋 Usage Workflow

### Patient information
- Open the accordion menu on the top to insert patient's data. Name, accession number and age are completely optional. Patient's sex is mandatory to obtain a correct Agatston score normalized to sex.

### Tab 0: Settings
- Install required dependencies (reportlab/matplotlib)
- Configure threshold (default: 130 HU for non-contrast CT)
- **Set company logo and description** for PDF reports
- View plugin information and version

### Tab 1: Setup View
1. Select CT volume (non-contrast, ECG-gated preferred). If a node is already active, the field will be automatically populated
2. (Optional) Apply optimal layout for better visualization
3. Tab 2 will enable automatically after volume selection

### Tab 2: Calcium Segmentation
1. Choose segmentation method:
   - **Point and Click** ⭐ *Recommended*: Click on calcium to automatically segment. After the selection enable the cursor with the button "Enable Point and click" and follow the instructions in the dialog box
   - **ROI Box**: Place a box around the aortic valve, use the cursor on the left to rotate the ROI according to the Aortic bulb. Please be careful not to include the coronary ostia → Apply threshold
   - **Paint Mode**: Manually paint calcium regions with threshold-guided brush. Click the erase button to delete any segmentation you don't need to include
2. The plugin automatically filters noise (lesions < 1.0 mm³)
3. Tab 3 will enable after segmentation

### Tab 3: Results & Analysis
1. Click **Calculate Agatston Score**
2. View quantitative results (Agatston Score, volume, equivalent mass, density, severity)
3. **Show 3D**: Visualize calcium with color-coded density levels
4. **Show Charts**: View histograms, percentiles, distributions per axial slice and bull's eye diagram
5. Generate all the data you would like to include in your PDF report
6. Tab 4 will enable after calculation

### Tab 4: Generate Report
1. Select output directory
2. Configure report options (screenshots, 3D view, charts)
3. Click **Generate PDF Report**
4. Final PDF report including your logo, patient info, results, clinical disclaimer, and references

## 🔬 Technical Details

### Agatston Score Calculation

```
Agatston Score = Σ (Lesion Area × Density Factor)
```

**Density Factors:** 1 (130-199 HU), 2 (200-299 HU), 3 (300-399 HU), 4 (≥400 HU)

**Noise Filtering:** Minimum lesion volume 1.5 mm³

### Severity Classification

**Women:** Severe ≥1377 AU | Moderate 400-1377 AU | Mild 100-400 AU  
**Men:** Severe ≥2062 AU | Moderate 1000-2062 AU | Mild 100-1000 AU

## 📚 References

1. **Agatston AS, et al.** Quantification of coronary artery calcium using ultrafast computed tomography. *J Am Coll Cardiol.* 1990;15(4):827-832.

2. **Clavel MA, et al.** Impact of Aortic Valve Calcification, as Measured by MDCT, on Survival in Patients With Aortic Stenosis. *JACC.* 2014;64(12):1202-1213.

3. **Pawade T, et al.** Computed Tomography Aortic Valve Calcium Scoring in Patients With Aortic Stenosis. *Circ Cardiovasc Imaging.* 2018;11(3):e007146.

4. **Tastet L, et al.** Sex-Related Differences in the Extent of Myocardial Fibrosis in Patients With Aortic Valve Stenosis. *JAHA.* 2020;9(8):e015344.

5. **ACC/AHA/AATS/ACR/ASA/SCA/SCAI/SIR/STS/SVM** Guidelines for the Management of Patients With Valvular Heart Disease. *Circulation.* 2021;143(5):e72-e227.

## 🤝 Contributing

Contributions welcome! Fork, create feature branch, commit, push, and open a Pull Request.

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Developed by**: Vittorio Censullo
- **Co-authored by**: **AITeRTC** - *Associazione Italiana Tecnici Esperti in TC*
- Built using **3D Slicer** open-source platform
- Based on **ACC/AHA** clinical guidelines

## 📞 Contact

**Vittorio Censullo**  
**AITeRTC** - Associazione Italiana Tecnici Esperti in TC  
GitHub: [@vcensullo](https://github.com/vcensullo)

Issues: [GitHub Issues](https://github.com/vcensullo/AoCalciumScore/issues)

If you use this plugin in your research please cite  

## 📖 Citation

```bibtex
@software{censullo2025aocalciumscore,
  author = {Censullo, Vittorio and AITeRTC},
  title = {AoCalciumScore: Aortic Valve Calcium Scoring Plugin for 3D Slicer},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/vcensullo/AoCalciumScore},
  institution = {AITeRTC - Associazione Italiana Tecnici di Radiologia e Radiodiagnostica e Tomografia Computerizzata ETS}
}
```

## 🔄 Version History

**Version 0.5** (October 2025) - Current
- Added Point and Click segmentation method
- Implemented settings persistence (QSettings)
- PDF report design with company branding
- Expanded scientific references
- Improved UI/UX throughout

---

**Version:** 0.5  
**Last Updated:** October 2025  
**Maintained by:** Vittorio Censullo & AITeRTC

*Vibe-coded with the aid of Claude*
