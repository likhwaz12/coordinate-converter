#!/usr/bin/env python3
"""
Coordinate Conversion Script
Converts coordinates from projected coordinate systems (like UTM) to latitude/longitude.
Reads from Excel files and outputs converted coordinates.
"""

import pandas as pd
import pyproj
from pyproj import Transformer
import argparse
import sys
import os
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Tuple
import xml.etree.ElementTree as ET
from xml.dom import minidom
import requests
import time
import zipfile
import math

class ElevationService:
    """Handles elevation data retrieval from various sources"""
    
    def __init__(self):
        self.cache = {}
        self.api_delay = 0.1  # Delay between API calls to be respectful
    
    def get_elevation_batch(self, coordinates: List[Tuple[float, float]], 
                          max_batch_size: int = 100) -> List[float]:
        """
        Get elevation data for multiple coordinates using Open Elevation API
        
        Args:
            coordinates: List of (latitude, longitude) tuples
            max_batch_size: Maximum number of coordinates per API call
            
        Returns:
            List of elevation values in meters
        """
        elevations = []
        
        for i in range(0, len(coordinates), max_batch_size):
            batch = coordinates[i:i + max_batch_size]
            batch_elevations = self._get_elevation_batch_api(batch)
            elevations.extend(batch_elevations)
            
            # Add delay between batches to be respectful to the API
            if i + max_batch_size < len(coordinates):
                time.sleep(self.api_delay)
        
        return elevations
    
    def _get_elevation_batch_api(self, coordinates: List[Tuple[float, float]]) -> List[float]:
        """Get elevation data from Open Elevation API"""
        try:
            # Prepare locations for API
            locations = []
            for lat, lon in coordinates:
                cache_key = f"{lat:.6f},{lon:.6f}"
                if cache_key in self.cache:
                    locations.append(None)  # Will be filled from cache
                else:
                    locations.append({"latitude": lat, "longitude": lon})
            
            # Get uncached elevations from API
            uncached_locations = [loc for loc in locations if loc is not None]
            api_elevations = []
            
            if uncached_locations:
                url = "https://api.open-elevation.com/api/v1/lookup"
                response = requests.post(url, json={"locations": uncached_locations}, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    api_elevations = [result["elevation"] for result in data["results"]]
                else:
                    # Fallback to zero elevation if API fails
                    print(f"Elevation API failed (status {response.status_code}), using ground level (0m)")
                    api_elevations = [0.0] * len(uncached_locations)
            
            # Combine cached and API results
            result_elevations = []
            api_index = 0
            
            for i, (lat, lon) in enumerate(coordinates):
                cache_key = f"{lat:.6f},{lon:.6f}"
                if cache_key in self.cache:
                    result_elevations.append(self.cache[cache_key])
                else:
                    elevation = api_elevations[api_index] if api_index < len(api_elevations) else 0.0
                    self.cache[cache_key] = elevation
                    result_elevations.append(elevation)
                    api_index += 1
            
            return result_elevations
            
        except Exception as e:
            print(f"Error getting elevation data: {e}")
            print("Using ground level (0m) for all turbines")
            return [0.0] * len(coordinates)
    
    def get_elevation_simple(self, lat: float, lon: float) -> float:
        """Get elevation for a single coordinate (uses batch method)"""
        elevations = self.get_elevation_batch([(lat, lon)])
        return elevations[0] if elevations else 0.0

@dataclass
class WindTurbineSpec:
    """Wind turbine 3D model specifications"""
    hub_height: float = 100.0  # meters - height to nacelle center
    rotor_diameter: float = 120.0  # meters - blade tip to tip
    tower_base_diameter: float = 6.0  # meters - tower base width
    tower_top_diameter: float = 3.5  # meters - tower top width
    nacelle_length: float = 15.0  # meters
    nacelle_width: float = 4.0  # meters
    nacelle_height: float = 4.5  # meters
    blade_length: float = 60.0  # meters - from hub to tip
    blade_chord: float = 3.0  # meters - blade width at base
    tower_color: str = "ffffffff"  # AABBGGRR format (white)
    nacelle_color: str = "ffcccccc"  # AABBGGRR format (light gray)
    blade_color: str = "ff999999"  # AABBGGRR format (gray)

class KMLGenerator:
    """Handles KML file generation for Google Earth visualization"""
    
    def __init__(self, turbine_spec: Optional[WindTurbineSpec] = None):
        self.turbine_spec = turbine_spec or WindTurbineSpec()
        self.kml_namespace = "http://www.opengis.net/kml/2.2"
        self.gx_namespace = "http://www.google.com/kml/ext/2.2"
    
    def create_kml_document(self, name: str = "Wind Farm", description: str = "3D Wind Turbines") -> ET.Element:
        """Create base KML document structure"""
        kml = ET.Element("kml")
        kml.set("xmlns", self.kml_namespace)
        kml.set("xmlns:gx", self.gx_namespace)
        
        document = ET.SubElement(kml, "Document")
        ET.SubElement(document, "name").text = name
        ET.SubElement(document, "description").text = description
        
        # Add styles for different turbine components
        self._add_turbine_styles(document)
        
        return kml
    
    def _add_turbine_styles(self, document: ET.Element):
        """Add KML styles for turbine components"""
        # Tower style
        tower_style = ET.SubElement(document, "Style")
        tower_style.set("id", "towerStyle")
        poly_style = ET.SubElement(tower_style, "PolyStyle")
        ET.SubElement(poly_style, "color").text = self.turbine_spec.tower_color
        ET.SubElement(poly_style, "fill").text = "1"
        line_style = ET.SubElement(tower_style, "LineStyle")
        ET.SubElement(line_style, "color").text = self.turbine_spec.tower_color
        ET.SubElement(line_style, "width").text = "2"
        
        # Nacelle style
        nacelle_style = ET.SubElement(document, "Style")
        nacelle_style.set("id", "nacelleStyle")
        poly_style = ET.SubElement(nacelle_style, "PolyStyle")
        ET.SubElement(poly_style, "color").text = self.turbine_spec.nacelle_color
        ET.SubElement(poly_style, "fill").text = "1"
        
        # Blade style
        blade_style = ET.SubElement(document, "Style")
        blade_style.set("id", "bladeStyle")
        poly_style = ET.SubElement(blade_style, "PolyStyle")
        ET.SubElement(poly_style, "color").text = self.turbine_spec.blade_color
        ET.SubElement(poly_style, "fill").text = "1"
    
    def create_tower_geometry(self, lat: float, lon: float, elevation: float = 0) -> List[List[Tuple[float, float, float]]]:
        """Generate improved tower geometry as a vertical tapered cylinder"""
        # Create tower as 12-sided polygon for good balance of quality and performance
        num_sides = 12
        angles = [i * (360 / num_sides) for i in range(num_sides)]
        
        base_radius = self.turbine_spec.tower_base_diameter / 2
        top_radius = self.turbine_spec.tower_top_diameter / 2
        hub_height = elevation + self.turbine_spec.hub_height
        
        # Calculate precise coordinate conversion factors
        lat_to_meter = 111320.0  # meters per degree latitude
        lon_to_meter = 111320.0 * math.cos(math.radians(lat))  # meters per degree longitude
        
        # Generate base and top circles in horizontal plane
        base_circle = []
        top_circle = []
        
        for angle in angles:
            angle_rad = math.radians(angle)
            
            # Base circle points (at ground elevation)
            base_dx = base_radius * math.cos(angle_rad) / lat_to_meter
            base_dy = base_radius * math.sin(angle_rad) / lon_to_meter
            base_circle.append((lon + base_dy, lat + base_dx, elevation))
            
            # Top circle points (at hub height)
            top_dx = top_radius * math.cos(angle_rad) / lat_to_meter
            top_dy = top_radius * math.sin(angle_rad) / lon_to_meter
            top_circle.append((lon + top_dy, lat + top_dx, hub_height))
        
        # Create tower wall polygons (vertical surfaces)
        tower_polygons = []
        
        for i in range(num_sides):
            next_i = (i + 1) % num_sides
            
            # Create trapezoidal panel between base and top (vertical wall)
            # Order vertices counterclockwise when viewed from outside
            wall_polygon = [
                base_circle[i],      # bottom-left
                top_circle[i],       # top-left
                top_circle[next_i],  # top-right
                base_circle[next_i], # bottom-right
                base_circle[i]       # close polygon
            ]
            tower_polygons.append(wall_polygon)
        
        return tower_polygons
    
    def create_nacelle_geometry(self, lat: float, lon: float, elevation: float = 0) -> List[List[Tuple[float, float, float]]]:
        """Generate nacelle positioned correctly - front third above tower, rest extending eastward"""
        hub_height = elevation + self.turbine_spec.hub_height
        length = self.turbine_spec.nacelle_length
        half_width = self.turbine_spec.nacelle_width / 2
        half_height = self.turbine_spec.nacelle_height / 2
        
        # Position nacelle so front third is above tower, rest extends eastward
        front_offset = length / 3  # Front third above tower
        back_extension = 2 * length / 3  # Back two-thirds extend eastward
        
        # Calculate coordinate conversion factors
        lat_to_meter = 111320.0
        lon_to_meter = 111320.0 * math.cos(math.radians(lat))
        
        # Nacelle extends eastward (positive longitude direction)
        # Front of nacelle (1/3 length) is centered over tower
        # Back of nacelle (2/3 length) extends east from tower center
        front_x = -front_offset / 2  # Start slightly west of tower center
        back_x = back_extension * 1.5   # Extend well east of tower center
        
        # Create 8 box vertices - positioned horizontally extending eastward
        vertices = []
        
        # Define the nacelle box corners
        for x in [front_x, back_x]:  # Front to back (west to east)
            for y in [-half_width, half_width]:  # Left to right (north to south)
                for z in [-half_height, half_height]:  # Bottom to top
                    lat_offset = x / lat_to_meter
                    lon_offset = y / lon_to_meter
                    vertices.append((lon + lon_offset, lat + lat_offset, hub_height + z))
        
        # Vertices order: [front-left-bottom, front-left-top, front-right-bottom, front-right-top,
        #                  back-left-bottom, back-left-top, back-right-bottom, back-right-top]
        box_faces = []
        
        # Front face (west face - over tower)
        front_face = [vertices[0], vertices[2], vertices[3], vertices[1], vertices[0]]
        box_faces.append(front_face)
        
        # Back face (east face - hanging in air)
        back_face = [vertices[4], vertices[5], vertices[7], vertices[6], vertices[4]]
        box_faces.append(back_face)
        
        # Bottom face
        bottom_face = [vertices[0], vertices[4], vertices[6], vertices[2], vertices[0]]
        box_faces.append(bottom_face)
        
        # Top face
        top_face = [vertices[1], vertices[3], vertices[7], vertices[5], vertices[1]]
        box_faces.append(top_face)
        
        # Left face (north face)
        left_face = [vertices[0], vertices[1], vertices[5], vertices[4], vertices[0]]
        box_faces.append(left_face)
        
        # Right face (south face)
        right_face = [vertices[2], vertices[6], vertices[7], vertices[3], vertices[2]]
        box_faces.append(right_face)
        
        return box_faces
    
    def create_blade_geometry(self, lat: float, lon: float, elevation: float = 0, 
                            blade_angle: float = 0) -> List[List[Tuple[float, float, float]]]:
        """Generate vertical blade geometry facing east (perpendicular to ground)"""
        hub_height = elevation + self.turbine_spec.hub_height
        blade_length = self.turbine_spec.blade_length
        blade_chord = self.turbine_spec.blade_chord
        
        # Calculate coordinate conversion factors
        lat_to_meter = 111320.0
        lon_to_meter = 111320.0 * math.cos(math.radians(lat))
        
        # Blades rotate in vertical plane, facing east
        # blade_angle rotates around the hub (0°, 120°, 240°)
        angle_rad = math.radians(blade_angle)
        
        # Create blade as vertical rectangular segments
        num_segments = 5
        blade_polygons = []
        
        for i in range(num_segments):
            # Distance ratios along blade (from hub to tip)
            start_ratio = i / num_segments
            end_ratio = (i + 1) / num_segments
            
            # Distances from hub center
            start_distance = start_ratio * blade_length
            end_distance = end_ratio * blade_length
            
            # Chord widths (tapered from hub to tip)
            start_chord = blade_chord * (1 - start_ratio * 0.7)  # Taper to 30% at tip
            end_chord = blade_chord * (1 - end_ratio * 0.7)
            
            # Calculate blade positions in vertical plane
            # Blades extend radially from hub in vertical plane
            # For vertical orientation: one component is height (Z), other is horizontal (Y or X)
            start_vertical = start_distance * math.cos(angle_rad)  # Vertical displacement from hub
            start_horizontal = start_distance * math.sin(angle_rad)  # Horizontal displacement in rotor plane
            end_vertical = end_distance * math.cos(angle_rad)
            end_horizontal = end_distance * math.sin(angle_rad)
            
            # Hub is at lat,lon,hub_height - blades extend in vertical plane
            # For east-facing turbine: blades sweep in north-south and up-down directions
            
            # Start section vertices (closer to hub)
            start_height = hub_height + start_vertical
            start_north = lat + start_horizontal / lat_to_meter
            
            # End section vertices (closer to tip)
            end_height = hub_height + end_vertical  
            end_north = lat + end_horizontal / lat_to_meter
            
            # Blade thickness extends eastward (into wind direction)
            start_thickness = start_chord / 2
            end_thickness = end_chord / 2
            
            # Create blade segment as rectangular panel in vertical plane
            # Each segment has 4 corners: front-hub, back-hub, back-tip, front-tip
            
            # Front edge (leading edge, facing wind from east)
            start_front = (lon + start_thickness / lon_to_meter, start_north, start_height)
            end_front = (lon + end_thickness / lon_to_meter, end_north, end_height)
            
            # Back edge (trailing edge)
            start_back = (lon - start_thickness / lon_to_meter, start_north, start_height)
            end_back = (lon - end_thickness / lon_to_meter, end_north, end_height)
            
            # Create front surface (wind-facing side)
            front_surface = [start_front, end_front, end_back, start_back, start_front]
            blade_polygons.append(front_surface)
            
            # Create back surface (downwind side)  
            back_surface = [start_back, end_back, end_front, start_front, start_back]
            blade_polygons.append(back_surface)
            
            # Create edge surfaces if not the last segment
            if i < num_segments - 1:
                # Leading edge surface
                next_start_ratio = (i + 1) / num_segments
                next_start_distance = next_start_ratio * blade_length
                next_start_y = next_start_distance * math.cos(angle_rad)
                next_start_z = next_start_distance * math.sin(angle_rad)
                next_start_height = hub_height + next_start_z
                next_start_north = lat + next_start_y / lat_to_meter
                next_start_thickness = blade_chord * (1 - next_start_ratio * 0.7) / 2
                
                next_start_front = (lon + next_start_thickness / lon_to_meter, next_start_north, next_start_height)
                next_start_back = (lon - next_start_thickness / lon_to_meter, next_start_north, next_start_height)
                
                # Leading edge panel
                leading_edge = [start_front, start_back, next_start_back, next_start_front, start_front]
                blade_polygons.append(leading_edge)
        
        return blade_polygons
    
    def add_turbine_placemark(self, parent: ET.Element, turbine_id: str, lat: float, lon: float, 
                            elevation: float = 0, description: str = "", use_3d_model: bool = True):
        """Add a complete wind turbine placemark to KML using COLLADA 3D model or basic geometry"""
        if use_3d_model:
            self._add_turbine_3d_model(parent, turbine_id, lat, lon, elevation, description)
        else:
            self._add_turbine_basic_geometry(parent, turbine_id, lat, lon, elevation, description)
    
    def _add_turbine_3d_model(self, parent: ET.Element, turbine_id: str, lat: float, lon: float, 
                            elevation: float = 0, description: str = ""):
        """Add wind turbine using COLLADA 3D model"""
        placemark = ET.SubElement(parent, "Placemark")
        ET.SubElement(placemark, "name").text = f"Turbine {turbine_id}"
        ET.SubElement(placemark, "description").text = description or f"Wind Turbine {turbine_id}\nHub Height: {self.turbine_spec.hub_height}m\nRotor Diameter: {self.turbine_spec.rotor_diameter}m"
        
        # Create Model element for COLLADA 3D model
        model = ET.SubElement(placemark, "Model")
        ET.SubElement(model, "altitudeMode").text = "absolute"
        
        # Position the model
        location = ET.SubElement(model, "Location")
        ET.SubElement(location, "longitude").text = str(lon)
        ET.SubElement(location, "latitude").text = str(lat)
        ET.SubElement(location, "altitude").text = str(elevation)
        
        # Orientation - face east (rotate 180° from current west-facing)
        orientation = ET.SubElement(model, "Orientation")
        ET.SubElement(orientation, "heading").text = "270"  # Face east (270° from north = 90° + 180°)
        ET.SubElement(orientation, "tilt").text = "0"      # Upright
        ET.SubElement(orientation, "roll").text = "0"      # No roll
        
        # Scale the model based on turbine specifications
        # Assume the DAE model is designed for 100m hub height and 120m rotor diameter
        scale_factor = self.turbine_spec.hub_height / 100.0
        scale = ET.SubElement(model, "Scale")
        ET.SubElement(scale, "x").text = str(scale_factor)
        ET.SubElement(scale, "y").text = str(scale_factor)
        ET.SubElement(scale, "z").text = str(scale_factor)
        
        # Link to the COLLADA model file
        link = ET.SubElement(model, "Link")
        ET.SubElement(link, "href").text = "models/turbine.dae"
        
    def _add_turbine_basic_geometry(self, parent: ET.Element, turbine_id: str, lat: float, lon: float, 
                                  elevation: float = 0, description: str = ""):
        """Add wind turbine using basic geometry (fallback method)"""
        # Create folder to organize turbine components
        folder = ET.SubElement(parent, "Folder")
        ET.SubElement(folder, "name").text = f"Turbine {turbine_id}"
        ET.SubElement(folder, "description").text = description or f"Wind Turbine {turbine_id}\nHub Height: {self.turbine_spec.hub_height}m\nRotor Diameter: {self.turbine_spec.rotor_diameter}m"
        
        # Add tower as separate placemark
        tower_placemark = ET.SubElement(folder, "Placemark")
        ET.SubElement(tower_placemark, "name").text = f"Tower {turbine_id}"
        self._add_tower_multigeometry(tower_placemark, lat, lon, elevation)
        
        # Add nacelle as separate placemark
        nacelle_placemark = ET.SubElement(folder, "Placemark")
        ET.SubElement(nacelle_placemark, "name").text = f"Nacelle {turbine_id}"
        self._add_nacelle_multigeometry(nacelle_placemark, lat, lon, elevation)
        
        # Add three blades as separate placemarks at 120-degree intervals
        for blade_num in range(3):
            blade_angle = blade_num * 120
            blade_placemark = ET.SubElement(folder, "Placemark")
            ET.SubElement(blade_placemark, "name").text = f"Blade {blade_num + 1} - Turbine {turbine_id}"
            self._add_blade_multigeometry(blade_placemark, lat, lon, elevation, blade_angle)
    
    def _add_tower_multigeometry(self, placemark: ET.Element, lat: float, lon: float, elevation: float):
        """Add improved tower geometry to placemark"""
        tower_polygons = self.create_tower_geometry(lat, lon, elevation)
        
        multigeometry = ET.SubElement(placemark, "MultiGeometry")
        
        for polygon_coords in tower_polygons:
            polygon = ET.SubElement(multigeometry, "Polygon")
            ET.SubElement(polygon, "altitudeMode").text = "absolute"
            ET.SubElement(polygon, "tessellate").text = "1"
            
            outer_boundary = ET.SubElement(polygon, "outerBoundaryIs")
            linear_ring = ET.SubElement(outer_boundary, "LinearRing")
            
            # Convert polygon coordinates to KML format
            coord_text = " ".join([f"{coord[0]},{coord[1]},{coord[2]}" for coord in polygon_coords])
            ET.SubElement(linear_ring, "coordinates").text = coord_text
            ET.SubElement(polygon, "styleUrl").text = "#towerStyle"
    
    def _add_nacelle_multigeometry(self, placemark: ET.Element, lat: float, lon: float, elevation: float):
        """Add improved nacelle geometry to placemark"""
        nacelle_faces = self.create_nacelle_geometry(lat, lon, elevation)
        
        multigeometry = ET.SubElement(placemark, "MultiGeometry")
        
        for face_coords in nacelle_faces:
            polygon = ET.SubElement(multigeometry, "Polygon")
            ET.SubElement(polygon, "altitudeMode").text = "absolute"
            ET.SubElement(polygon, "tessellate").text = "1"
            
            outer_boundary = ET.SubElement(polygon, "outerBoundaryIs")
            linear_ring = ET.SubElement(outer_boundary, "LinearRing")
            
            # Convert face coordinates to KML format
            coord_text = " ".join([f"{coord[0]},{coord[1]},{coord[2]}" for coord in face_coords])
            ET.SubElement(linear_ring, "coordinates").text = coord_text
            ET.SubElement(polygon, "styleUrl").text = "#nacelleStyle"
    
    def _add_blade_multigeometry(self, placemark: ET.Element, lat: float, lon: float, 
                               elevation: float, blade_angle: float):
        """Add improved blade geometry to placemark"""
        blade_polygons = self.create_blade_geometry(lat, lon, elevation, blade_angle)
        
        multigeometry = ET.SubElement(placemark, "MultiGeometry")
        
        for polygon_coords in blade_polygons:
            polygon = ET.SubElement(multigeometry, "Polygon")
            ET.SubElement(polygon, "altitudeMode").text = "absolute"
            ET.SubElement(polygon, "tessellate").text = "1"
            
            outer_boundary = ET.SubElement(polygon, "outerBoundaryIs")
            linear_ring = ET.SubElement(outer_boundary, "LinearRing")
            
            # Convert polygon coordinates to KML format
            coord_text = " ".join([f"{coord[0]},{coord[1]},{coord[2]}" for coord in polygon_coords])
            ET.SubElement(linear_ring, "coordinates").text = coord_text
            ET.SubElement(polygon, "styleUrl").text = "#bladeStyle"
    
    def save_kml(self, kml_element: ET.Element, filepath: Path):
        """Save KML to file with proper formatting"""
        rough_string = ET.tostring(kml_element, 'unicode')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")
        
        # Remove empty lines
        pretty_xml = '\n'.join(line for line in pretty_xml.split('\n') if line.strip())
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(pretty_xml)
    
    def add_wind_farm_camera_view(self, document: ET.Element, coordinates: List[Tuple[float, float, float]]):
        """Add LookAt element to automatically zoom to wind farm area when KMZ opens"""
        if not coordinates:
            return
        
        # Calculate bounding box of all turbine coordinates
        lats = [coord[0] for coord in coordinates]
        lons = [coord[1] for coord in coordinates]
        elevs = [coord[2] for coord in coordinates]
        
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
        min_elev, max_elev = min(elevs), max(elevs)
        
        # Calculate center point and size of wind farm
        center_lat = (min_lat + max_lat) / 2
        center_lon = (min_lon + max_lon) / 2
        center_elev = (min_elev + max_elev) / 2
        
        # Calculate span in degrees
        lat_span = max_lat - min_lat
        lon_span = max_lon - min_lon
        
        # Calculate appropriate viewing distance
        # Convert degree spans to approximate kilometers
        lat_km = lat_span * 111.32  # 1 degree lat ≈ 111.32 km
        lon_km = lon_span * 111.32 * math.cos(math.radians(center_lat))  # longitude distance varies with latitude
        
        # Use the larger span to determine range, add buffer for good framing
        max_span_km = max(lat_km, lon_km, 1.0)  # Minimum 1km span
        range_meters = max_span_km * 1000 * 2.5  # 2.5x buffer for good framing
        
        # Add hub height for better 3D view perspective
        altitude = center_elev + self.turbine_spec.hub_height
        
        # Create LookAt element for automatic camera positioning
        lookat = ET.SubElement(document, "LookAt")
        ET.SubElement(lookat, "longitude").text = str(center_lon)
        ET.SubElement(lookat, "latitude").text = str(center_lat)
        ET.SubElement(lookat, "altitude").text = str(altitude)
        ET.SubElement(lookat, "heading").text = "0"      # North orientation
        ET.SubElement(lookat, "tilt").text = "45"        # 45-degree tilt for good 3D perspective
        ET.SubElement(lookat, "range").text = str(range_meters)  # Distance from viewpoint to LookAt position
        ET.SubElement(lookat, "altitudeMode").text = "absolute"
        
        print(f"📍 Auto-zoom configured:")
        print(f"   Center: {center_lat:.6f}, {center_lon:.6f}")
        print(f"   Wind farm span: {lat_km:.2f}km × {lon_km:.2f}km")
        print(f"   Viewing range: {range_meters/1000:.2f}km")

    def save_kmz(self, kml_element: ET.Element, filepath: Path, include_3d_models: bool = True):
        """Save KML as compressed KMZ file with 3D models"""
        # Create temporary KML file
        kml_temp_path = filepath.with_suffix('.kml')
        self.save_kml(kml_element, kml_temp_path)
        
        # Create KMZ (ZIP archive containing KML and models)
        with zipfile.ZipFile(filepath, 'w', zipfile.ZIP_DEFLATED) as kmz:
            # Add the KML file as doc.kml (required name for KMZ)
            kmz.write(kml_temp_path, 'doc.kml')
            
            # Add 3D model files if requested and available
            if include_3d_models:
                turbine_dae_path = Path(__file__).parent / 'turbine.dae'
                if turbine_dae_path.exists():
                    # Add the COLLADA model to models subdirectory in KMZ
                    kmz.write(turbine_dae_path, 'models/turbine.dae')
                    print(f"✅ Included 3D model: turbine.dae")
                else:
                    print(f"⚠️  Warning: turbine.dae not found at {turbine_dae_path}")
                    print(f"   Falling back to basic geometry")
        
        # Clean up temporary KML file
        kml_temp_path.unlink()
        
        print(f"📦 KMZ file created: {filepath}")

def detect_coordinate_system(easting, northing):
    """
    Attempt to detect the coordinate system based on coordinate values.
    Returns the most likely EPSG code.
    """
    # Check for UTM zones based on coordinate ranges
    if 166000 <= easting <= 834000:  # Typical UTM easting range
        if 0 <= northing <= 10000000:  # Northern hemisphere
            # Estimate UTM zone from easting
            zone = int((easting - 166000) / 1000000) + 1
            if 1 <= zone <= 60:
                # Assume WGS84 UTM for now - can be made more sophisticated
                epsg_code = 32600 + zone  # WGS84 UTM North
                return f"EPSG:{epsg_code}", f"WGS84 UTM Zone {zone}N"
        elif northing > 10000000:  # Southern hemisphere
            zone = int((easting - 166000) / 1000000) + 1
            if 1 <= zone <= 60:
                epsg_code = 32700 + zone  # WGS84 UTM South
                return f"EPSG:{epsg_code}", f"WGS84 UTM Zone {zone}S"
    
    # Default fallback - assume UTM Zone 36S based on your example
    return "EPSG:32736", "WGS84 UTM Zone 36S"

def convert_coordinates(input_file, output_file=None, source_crs="EPSG:32736", target_crs="EPSG:4326"):
    """
    Convert coordinates from source CRS to target CRS.
    
    Args:
        input_file: Path to input Excel file (relative to Inputs folder or absolute path)
        output_file: Path to output file (optional, will be saved in Outputs folder)
        source_crs: Source coordinate reference system (default: EPSG:32736 - WGS84 UTM Zone 36S)
        target_crs: Target coordinate reference system (default: WGS84 lat/lon)
    """
    
    try:
        # Create Inputs and Outputs directories if they don't exist
        inputs_dir = Path("Inputs")
        outputs_dir = Path("Outputs")
        inputs_dir.mkdir(exist_ok=True)
        outputs_dir.mkdir(exist_ok=True)
        
        # Construct full input path
        input_path = Path(input_file)
        if not input_path.is_absolute():
            input_path = inputs_dir / input_file
        
        # Check if input file exists
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")
        
        # Read the Excel file
        print(f"Reading Excel file: {input_path}")
        df = pd.read_excel(input_path)
        
        # Display the first few rows and column names
        print("\nInput data structure:")
        print(df.head())
        print(f"\nColumns: {list(df.columns)}")
        
        # Try to identify coordinate columns
        coord_columns = {}
        for col in df.columns:
            col_lower = col.lower()
            if any(word in col_lower for word in ['east', 'x', 'utm_x']):
                coord_columns['easting'] = col
            elif any(word in col_lower for word in ['north', 'south', 'y', 'utm_y']):
                coord_columns['northing'] = col
        
        # If standard names not found, assume common patterns
        if not coord_columns:
            possible_cols = list(df.columns)
            if len(possible_cols) >= 3:  # Assuming fid, easting, northing
                coord_columns['easting'] = possible_cols[1]
                coord_columns['northing'] = possible_cols[2]
            else:
                print("Error: Could not identify coordinate columns")
                return False
        
        print(f"Using columns - Easting: '{coord_columns['easting']}', Northing: '{coord_columns['northing']}'")
        
        # Get coordinate values
        easting = df[coord_columns['easting']].values
        northing = df[coord_columns['northing']].values
        
        # Auto-detect coordinate system if "auto" is specified, otherwise use provided CRS
        if source_crs == "auto":
            source_crs, crs_description = detect_coordinate_system(easting[0], northing[0])
            print(f"Auto-detected coordinate system: {crs_description} ({source_crs})")
        else:
            # Determine description for common EPSG codes
            crs_descriptions = {
                "EPSG:32736": "WGS84 UTM Zone 36S",
                "EPSG:32636": "WGS84 UTM Zone 36N",
                "EPSG:32737": "WGS84 UTM Zone 37S",
                "EPSG:32637": "WGS84 UTM Zone 37N",
                "EPSG:4326": "WGS84 Latitude/Longitude"
            }
            crs_description = crs_descriptions.get(source_crs, source_crs)
            print(f"Using specified coordinate system: {crs_description} ({source_crs})")
        
        # Create transformer
        print(f"Converting from {source_crs} to {target_crs}")
        transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
        
        # Convert coordinates
        print("Converting coordinates...")
        longitude, latitude = transformer.transform(easting, northing)
        
        # Create output dataframe
        output_df = df.copy()
        output_df['Latitude'] = latitude
        output_df['Longitude'] = longitude
        
        # Reorder columns to put lat/lon first after any ID column
        cols = list(output_df.columns)
        if 'fid' in cols or 'id' in [c.lower() for c in cols]:
            id_col = next((c for c in cols if c.lower() in ['fid', 'id']), cols[0])
            new_cols = [id_col, 'Latitude', 'Longitude'] + [c for c in cols if c not in [id_col, 'Latitude', 'Longitude']]
        else:
            new_cols = ['Latitude', 'Longitude'] + [c for c in cols if c not in ['Latitude', 'Longitude']]
        
        output_df = output_df[new_cols]
        
        # Determine output file name and path
        if output_file is None:
            output_filename = f"{input_path.stem}_converted.xlsx"
            output_path = outputs_dir / output_filename
        else:
            output_path = Path(output_file)
            if not output_path.is_absolute():
                output_path = outputs_dir / output_file
        
        # Save to Excel
        print(f"Saving converted coordinates to: {output_path}")
        output_df.to_excel(output_path, index=False)
        
        # Also save as CSV for convenience
        csv_file = output_path.with_suffix('.csv')
        output_df.to_csv(csv_file, index=False)
        print(f"Also saved as CSV: {csv_file}")
        
        # Display summary
        print(f"\nConversion Summary:")
        print(f"Input file: {input_path}")
        print(f"Output file: {output_path}")
        print(f"Records processed: {len(output_df)}")
        print(f"Source CRS: {source_crs}")
        print(f"Target CRS: {target_crs}")
        
        print(f"\nFirst few converted coordinates:")
        print(output_df[['Latitude', 'Longitude']].head())
        
        return True
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return False

def generate_google_earth_kml(csv_file_path, output_dir=None, turbine_spec=None, use_elevation=True, save_kmz=True, use_3d_models=True):
    """
    Generate Google Earth KML file from converted coordinates CSV with enhanced 3D models
    
    Args:
        csv_file_path: Path to CSV file with Latitude and Longitude columns
        output_dir: Output directory (defaults to Outputs folder)
        turbine_spec: WindTurbineSpec object for customization
        use_elevation: Whether to fetch real elevation data (default: True)
        save_kmz: Whether to also save as compressed KMZ (default: True)
        use_3d_models: Whether to use COLLADA 3D models instead of basic geometry (default: True)
    """
    try:
        # Setup paths
        csv_path = Path(csv_file_path)
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        
        if output_dir is None:
            output_dir = Path("Outputs")
        else:
            output_dir = Path(output_dir)
        
        output_dir.mkdir(exist_ok=True)
        
        # Read coordinate data
        print(f"Reading coordinate data from: {csv_path}")
        df = pd.read_csv(csv_path)
        
        # Validate required columns
        required_cols = ['Latitude', 'Longitude']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Get elevation data if requested
        elevations = []
        if use_elevation:
            print("🌍 Fetching elevation data for terrain-aware placement...")
            elevation_service = ElevationService()
            coordinates = [(row['Latitude'], row['Longitude']) for _, row in df.iterrows()]
            elevations = elevation_service.get_elevation_batch(coordinates, max_batch_size=50)
            print(f"✅ Retrieved elevation data for {len(elevations)} locations")
        else:
            elevations = [0.0] * len(df)
        
        # Initialize KML generator
        kml_gen = KMLGenerator(turbine_spec)
        
        # Create KML document
        kml_root = kml_gen.create_kml_document(
            name=f"Wind Farm - {csv_path.stem}",
            description=f"Enhanced 3D Wind Turbine visualization generated from {csv_path.name}"
        )
        
        document = kml_root.find(".//Document")
        
        # Collect turbine coordinates for camera view calculation
        turbine_coordinates = []
        for idx, row in df.iterrows():
            lat = row['Latitude']
            lon = row['Longitude']
            elevation = elevations[idx] if idx < len(elevations) else 0.0
            turbine_coordinates.append((lat, lon, elevation))
        
        # Add automatic camera view to zoom to wind farm area BEFORE adding placemarks
        print(f"📹 Setting up automatic zoom to wind farm...")
        kml_gen.add_wind_farm_camera_view(document, turbine_coordinates)
        
        # Add turbines to KML with elevation data
        print(f"🏗️  Generating enhanced 3D turbines for {len(df)} locations...")
        for idx, row in df.iterrows():
            lat = row['Latitude']
            lon = row['Longitude']
            elevation = elevations[idx] if idx < len(elevations) else 0.0
            
            # Use fid or index as turbine ID
            turbine_id = row.get('fid', idx + 1)
            
            # Create enhanced description with available data
            description_parts = [f"Turbine ID: {turbine_id}"]
            if 'Easting' in df.columns and 'Southing' in df.columns:
                description_parts.append(f"Original Coordinates: {row['Easting']:.1f}, {row['Southing']:.1f}")
            description_parts.append(f"Lat/Lon: {lat:.6f}, {lon:.6f}")
            description_parts.append(f"Ground Elevation: {elevation:.1f}m")
            description_parts.append(f"Hub Height: {kml_gen.turbine_spec.hub_height}m")
            description_parts.append(f"Total Height: {elevation + kml_gen.turbine_spec.hub_height:.1f}m")
            description = "\n".join(description_parts)
            
            # Add turbine to KML with real elevation
            kml_gen.add_turbine_placemark(
                document, 
                str(turbine_id), 
                lat, 
                lon, 
                elevation=elevation,
                description=description,
                use_3d_model=use_3d_models
            )
        
        # Save KML file
        kml_filename = f"{csv_path.stem}_wind_farm_3d.kml"
        kml_path = output_dir / kml_filename
        
        print(f"💾 Saving enhanced KML file: {kml_path}")
        kml_gen.save_kml(kml_root, kml_path)
        
        # Save KMZ file if requested
        kmz_path = None
        if save_kmz:
            kmz_filename = f"{csv_path.stem}_wind_farm_3d.kmz"
            kmz_path = output_dir / kmz_filename
            print(f"📦 Saving compressed KMZ file: {kmz_path}")
            kml_gen.save_kmz(kml_root, kmz_path, include_3d_models=use_3d_models)
        
        # Display summary
        print(f"\n🎯 Enhanced Google Earth Generation Summary:")
        print(f"Input CSV: {csv_path}")
        print(f"Output KML: {kml_path}")
        if kmz_path:
            print(f"Output KMZ: {kmz_path}")
        print(f"Turbines generated: {len(df)}")
        print(f"Hub height: {kml_gen.turbine_spec.hub_height}m")
        print(f"Rotor diameter: {kml_gen.turbine_spec.rotor_diameter}m")
        print(f"Elevation data: {'✅ Real terrain data' if use_elevation else '❌ Ground level (0m)'}")
        if use_3d_models:
            print(f"3D Model quality: ✅ Professional COLLADA 3D models:")
            print(f"  • High-detail wind turbine geometry")
            print(f"  • Realistic proportions and textures")
            print(f"  • Optimized for Google Earth display")
        else:
            print(f"3D Model quality: ✅ Basic geometry:")
            print(f"  • 12-sided tapered tower cylinders")
            print(f"  • Realistic 6-faced nacelle boxes") 
            print(f"  • Vertical blade surfaces")
        print(f"\n✅ Enhanced KML/KMZ files ready for Google Earth!")
        print(f"📋 Instructions:")
        print(f"   1. Open Google Earth Pro")
        print(f"   2. File → Open → {kmz_path if kmz_path else kml_path}")
        print(f"   3. Navigate to turbine locations to view detailed 3D models")
        print(f"   4. Use 3D view and tilt for best visualization")
        
        return True
        
    except Exception as e:
        print(f"Error generating enhanced KML: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Convert coordinates from projected systems to latitude/longitude",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python coordinate_converter.py input.xlsx
  python coordinate_converter.py input.xlsx -o output.xlsx
  python coordinate_converter.py input.xlsx -s "EPSG:32636"
  python coordinate_converter.py input.xlsx -s "auto"
  python coordinate_converter.py input.xlsx -s "EPSG:32736" -t "EPSG:4326"
  
Note: Files are read from 'Inputs' folder and saved to 'Outputs' folder automatically.
These folders will be created if they don't exist.
Default source coordinate system is WGS84 UTM Zone 36S (EPSG:32736).
  
Common EPSG codes:
  EPSG:4326  - WGS84 Latitude/Longitude
  EPSG:32736 - WGS84 UTM Zone 36S (default)
  EPSG:32636 - WGS84 UTM Zone 36N
  EPSG:32737 - WGS84 UTM Zone 37S
  EPSG:32637 - WGS84 UTM Zone 37N
  EPSG:3857  - Web Mercator
        """
    )
    
    parser.add_argument('input_file', help='Input Excel file path (will look in Inputs folder if not absolute path)')
    parser.add_argument('-o', '--output', help='Output file path (will save to Outputs folder if not absolute path)')
    parser.add_argument('-s', '--source-crs', default='EPSG:32736',
                       help='Source coordinate reference system (default: EPSG:32736 - WGS84 UTM Zone 36S, use "auto" for auto-detection)')
    parser.add_argument('-t', '--target-crs', default='EPSG:4326', 
                       help='Target coordinate reference system (default: EPSG:4326 - WGS84 lat/lon)')
    parser.add_argument('--generate-kml', action='store_true',
                       help='Generate Google Earth KML file with 3D wind turbines')
    parser.add_argument('--hub-height', type=float, default=100.0,
                       help='Wind turbine hub height in meters (default: 100.0)')
    parser.add_argument('--rotor-diameter', type=float, default=120.0,
                       help='Wind turbine rotor diameter in meters (default: 120.0)')
    parser.add_argument('--no-elevation', action='store_true',
                       help='Skip elevation data fetching (use ground level)')
    parser.add_argument('--no-kmz', action='store_true',
                       help='Skip KMZ file generation (KML only)')
    parser.add_argument('--basic-geometry', action='store_true',
                       help='Use basic geometry instead of COLLADA 3D models')
    
    args = parser.parse_args()
    
    # The script will automatically look in Inputs folder and save to Outputs folder
    # No need to check if input file exists here since convert_coordinates handles it
    
    # Perform conversion
    success = convert_coordinates(
        input_file=args.input_file,
        output_file=args.output,
        source_crs=args.source_crs,
        target_crs=args.target_crs
    )
    
    if success:
        print("\n✅ Conversion completed successfully!")
        
        # Generate KML if requested
        if args.generate_kml:
            print("\n🌍 Generating Google Earth KML file...")
            
            # Create custom turbine specifications
            custom_turbine_spec = WindTurbineSpec(
                hub_height=args.hub_height,
                rotor_diameter=args.rotor_diameter,
                blade_length=args.rotor_diameter / 2  # blade length is half the rotor diameter
            )
            
            # Determine CSV file path from conversion
            inputs_dir = Path("Inputs")
            outputs_dir = Path("Outputs")
            
            input_path = Path(args.input_file)
            if not input_path.is_absolute():
                input_path = inputs_dir / args.input_file
            
            # Find the corresponding CSV file
            if args.output:
                csv_filename = Path(args.output).with_suffix('.csv').name
            else:
                csv_filename = f"{input_path.stem}_converted.csv"
            
            csv_path = outputs_dir / csv_filename
            
            # Generate enhanced KML with COLLADA 3D models
            kml_success = generate_google_earth_kml(
                csv_path,
                output_dir=outputs_dir,
                turbine_spec=custom_turbine_spec,
                use_elevation=not args.no_elevation,
                save_kmz=not args.no_kmz,
                use_3d_models=not args.basic_geometry
            )
            
            if not kml_success:
                print("\n❌ KML generation failed!")
                sys.exit(1)
    else:
        print("\n❌ Conversion failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()
