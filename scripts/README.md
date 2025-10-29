# TeloX Visualization Scripts

Standalone visualization tools for TeloX telomere annotation output.

## 🖥️ Web-Based Visualizer (JavaScript)

### **telomere_visualizer.html**
Interactive web-based visualization with drag-and-drop file upload.

**Usage:**
```bash
# Open in web browser
open telomere_visualizer.html
# or
firefox telomere_visualizer.html
```

**Features:**
- ✅ **No dependencies** - works in any modern browser
- ✅ **Drag & drop** file upload
- ✅ **Interactive filtering** by chromosome, motif, strand
- ✅ **Real-time statistics** 
- ✅ **Hover details** on all data points
- ✅ **Color-coded** by strand (red=forward, blue=reverse)

## 🐍 Python Visualizer

### **visualize_telomeres.py**
Comprehensive Python-based visualization with multiple output formats.

**Basic Usage:**
```bash
# Text summary + ASCII plots (always works)
python3 visualize_telomeres.py telox_genome_telomeres.bed

# With matplotlib plots (requires matplotlib)
python3 visualize_telomeres.py telox_genome_telomeres.bed

# Interactive HTML (requires plotly)
python3 visualize_telomeres.py telox_genome_telomeres.bed --interactive

# TSV input format
python3 visualize_telomeres.py telox_genome_annotations.tsv --format tsv
```

**Advanced Options:**
```bash
# Text-only mode (no dependencies)
python3 visualize_telomeres.py telox_genome_telomeres.bed --text-only

# Custom output prefix
python3 visualize_telomeres.py telox_genome_telomeres.bed --output-prefix my_genome

# Interactive mode
python3 visualize_telomeres.py telox_genome_telomeres.bed --interactive --output-prefix my_genome
```

## 📊 Output Files

### **Python Script Generates:**
- `telomere_ideogram.png` - Chromosome ideogram with telomere positions
- `telomere_summary.png` - 4-panel summary statistics
- `telomere_interactive.html` - Interactive Plotly visualization (with --interactive)

### **JavaScript Tool Generates:**
- Interactive plots directly in the browser
- No files saved (purely web-based)

## 🔧 Dependencies

### **JavaScript Tool:**
- **None** - works in any modern web browser
- Uses Plotly.js from CDN

### **Python Tool:**
- **Minimum**: Python 3.6+ (text mode only)
- **Recommended**: `pip install matplotlib numpy plotly`
- **Alternative**: `conda install matplotlib numpy plotly`

## 📋 Input Formats

### **BED Format (from TeloX --annotate):**
```
track name="TeloX_Telomeres" description="Telomere motifs identified by TeloX" itemRgb="On"
CM060242.1	169	174	GGTGT	+	255,0,0
CM060242.1	1250	1256	TTAGGG	-	0,0,255
```

### **TSV Format (from TeloX annotations.tsv):**
```
chromosome	start	end	motif	strand	score	region_type	length
CM060242.1	169	174	GGTGT	+	0.850	telomere	5
```

## 🎯 Recommended Workflow

### **Quick Start:**
```bash
# 1. Run TeloX
./target/release/telox genome.fasta --annotate

# 2. Visualize (choose one):
# Option A: Web browser (easiest)
open scripts/telomere_visualizer.html

# Option B: Python (most features)
python3 scripts/visualize_telomeres.py telox_genome_telomeres.bed
```

### **Publication Quality:**
```bash
# Generate high-resolution plots
python3 scripts/visualize_telomeres.py telox_genome_telomeres.bed --output-prefix publication

# Files generated:
# - publication_ideogram.png (300 DPI)
# - publication_summary.png (300 DPI)
```

### **Interactive Analysis:**
```bash
# Generate interactive HTML
python3 scripts/visualize_telomeres.py telox_genome_telomeres.bed --interactive

# Open telomere_interactive.html in browser for exploration
```

## 🎨 Visualization Features

### **Chromosome Ideogram:**
- Horizontal bars representing chromosomes
- Colored regions showing telomere positions
- Red = forward strand, Blue = reverse strand
- Telomere count annotations

### **Summary Statistics:**
- Chromosome distribution bar chart
- Top motifs frequency chart
- Strand bias pie chart
- Quality score histogram

### **Interactive Features:**
- Zoom and pan on all plots
- Hover for detailed information
- Filter by chromosome, motif, or strand
- Real-time statistics updates

## 🐛 Troubleshooting

### **No plots generated:**
```bash
# Check dependencies
python3 -c "import matplotlib; print('matplotlib OK')"
python3 -c "import plotly; print('plotly OK')"

# Use text-only mode
python3 visualize_telomeres.py file.bed --text-only
```

### **Empty visualization:**
```bash
# Check file format
head -5 telox_genome_telomeres.bed

# Try lower annotation threshold in TeloX
./target/release/telox genome.fasta --annotate --annotation-threshold 0.1
```

### **Browser issues:**
- Use modern browser (Chrome, Firefox, Safari, Edge)
- Enable JavaScript
- Check browser console for errors

## 💡 Tips

1. **Start with JavaScript tool** for quick exploration
2. **Use Python tool** for publication-quality figures
3. **Try text-only mode** if having dependency issues
4. **Adjust TeloX thresholds** if no annotations found
5. **Use interactive mode** for detailed analysis

The visualization tools are completely independent of TeloX and can be used with any compatible BED or TSV file!
