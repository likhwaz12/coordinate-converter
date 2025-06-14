# Coordinate Conversion Tool

This Python script converts coordinates from projected coordinate systems (like UTM) to latitude and longitude format. The script automatically organizes files using `Inputs` and `Outputs` folders.

## Folder Structure

The script expects and creates the following folder structure:
```
your_project/
├── coordinate_converter.py
├── requirements.txt
├── Inputs/
│   └── your_excel_files.xlsx
└── Outputs/
    ├── converted_files.xlsx
    └── converted_files.csv
```

## Setup

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Create folder structure:**
   The script will automatically create `Inputs` and `Outputs` folders if they don't exist.

3. **Place your Excel files in the `Inputs` folder**

## Usage

### Basic Usage (Default: WGS84 UTM Zone 36S)
Place your Excel file in the `Inputs` folder, then run:

```bash
python coordinate_converter.py your_file.xlsx
```

This will assume your coordinates are in **WGS84 UTM Zone 36S** (EPSG:32736) by default.

### Auto-detect Coordinate System
To let the script automatically detect the coordinate system:

```bash
python coordinate_converter.py your_file.xlsx -s "auto"
```

### Specify Different Source Coordinate System
```bash
python coordinate_converter.py your_file.xlsx -s "EPSG:32636"  # UTM Zone 36N
python coordinate_converter.py your_file.xlsx -s "EPSG:32737"  # UTM Zone 37S
```

### Specify Output File Name
```bash
python coordinate_converter.py input.xlsx -o custom_output_name.xlsx
```

### Complete Workflow Example
```bash
# 1. Place 50Nordexlayout.xlsx in the Inputs folder
# 2. Run the conversion (uses UTM Zone 36S by default)
python coordinate_converter.py 50Nordexlayout.xlsx
# 3. Find results in Outputs folder
```

**Note:** The script now defaults to **WGS84 UTM Zone 36S** (EPSG:32736). All input files are read from the `Inputs` folder, and all output files are saved to the `Outputs` folder automatically.

## Supported Input Formats

- **Excel files** (.xlsx, .xls)
- **Expected columns**: The script looks for columns containing coordinates:
  - Easting/X coordinates (columns with names like "Easting", "X", "UTM_X")
  - Northing/Y coordinates (columns with names like "Northing", "Southing", "Y", "UTM_Y")
- **Your file format**: Based on your Excel file, it expects columns like:
  - Column 1: ID (fid)
  - Column 2: Easting coordinates
  - Column 3: Northing coordinates

## Output

The script generates two output files in the `Outputs` folder:
1. **Excel file** (.xlsx) with original data plus Latitude and Longitude columns
2. **CSV file** (.csv) with the same data for easy import into other applications

Output files are automatically named as `{input_filename}_converted.xlsx` unless you specify a custom name.

## Common Coordinate Systems (EPSG Codes)

- `EPSG:4326` - WGS84 Latitude/Longitude (default output)
- `EPSG:32736` - **WGS84 UTM Zone 36S (default input)**
- `EPSG:32636` - WGS84 UTM Zone 36N
- `EPSG:32737` - WGS84 UTM Zone 37S  
- `EPSG:32637` - WGS84 UTM Zone 37N
- `EPSG:3857` - Web Mercator (used by Google Maps)

## Auto-Detection Feature

The script attempts to automatically detect UTM zones based on coordinate ranges:
- Analyzes Easting and Northing values
- Estimates the UTM zone
- Assumes WGS84 datum by default

## Example Workflow

1. **Setup your project:**
   ```bash
   mkdir coordinate_project
   cd coordinate_project
   # Place coordinate_converter.py and requirements.txt here
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Place your file:**
   ```bash
   # The script will create Inputs folder automatically, or create it manually:
   mkdir Inputs
   # Copy your Excel file to the Inputs folder
   cp /path/to/your/50Nordexlayout.xlsx Inputs/
   ```

4. **Run conversion:**
   ```bash
   # Default - assumes WGS84 UTM Zone 36S
   python coordinate_converter.py 50Nordexlayout.xlsx
   
   # Or with auto-detection
   python coordinate_converter.py 50Nordexlayout.xlsx -s "auto"
   ```

5. **Find results:**
   ```bash
   ls Outputs/
   # You'll see: 50Nordexlayout_converted.xlsx and 50Nordexlayout_converted.csv
   ```

## Troubleshooting

1. **"Input file not found"**
   - Make sure your Excel file is in the `Inputs` folder
   - Check the filename spelling and extension

2. **"Could not identify coordinate columns"**
   - Check that your Excel file has recognizable column names
   - Ensure columns contain numeric coordinate data

3. **Conversion errors**
   - Verify the source coordinate system is correct
   - Check that coordinates are in the expected format/range

4. **Permission errors**
   - Ensure you have write permissions for the `Outputs` folder
   - Close any open Excel files that might be locked

## Command Line Help

For complete usage information:
```bash
python coordinate_converter.py --help
```