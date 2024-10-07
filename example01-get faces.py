from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopoDS import topods
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Display.SimpleGui import init_display
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.AIS import AIS_Shape
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Lin
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeEdge

def read_step_file(file_path):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(file_path)
    if status != 1:
        print("Error: Cannot read STEP file")
        return None
    step_reader.TransferRoots()
    shape = step_reader.Shape()
    return shape

def display_origin(display):
    origin = gp_Pnt(0, 0, 0)
    vertex_O = BRepBuilderAPI_MakeVertex(origin).Vertex()
    ais_vertex_O = AIS_Shape(vertex_O)
    display.Context.Display(ais_vertex_O, False)
    display.DisplayMessage(origin, "O")

def display_axis_lines(display):
    origin = gp_Pnt(0, 0, 0)
    
    x_dir_line = gp_Lin(origin, gp_Dir(1, 0, 0))
    y_dir_line = gp_Lin(origin, gp_Dir(0, 1, 0))
    z_dir_line = gp_Lin(origin, gp_Dir(0, 0, 1))
    
    x_dir_edge = BRepBuilderAPI_MakeEdge(x_dir_line, -1000, 1000).Edge()
    y_dir_edge = BRepBuilderAPI_MakeEdge(y_dir_line, -1000, 1000).Edge()
    z_dir_edge = BRepBuilderAPI_MakeEdge(z_dir_line, -1000, 1000).Edge()
    
    ais_x_dir = AIS_Shape(x_dir_edge)
    ais_y_dir = AIS_Shape(y_dir_edge)
    ais_z_dir = AIS_Shape(z_dir_edge)
    
    yellow = Quantity_Color(1.0, 1.0, 0.0, Quantity_TOC_RGB)
    
    ais_x_dir.SetColor(yellow)
    ais_y_dir.SetColor(yellow)
    ais_z_dir.SetColor(yellow)
    
    display.Context.Display(ais_x_dir, False)
    display.Context.Display(ais_y_dir, False)
    display.Context.Display(ais_z_dir, False)
    
    display.DisplayMessage(gp_Pnt(1000, 0, 0), "x-dir")
    display.DisplayMessage(gp_Pnt(0, 1000, 0), "y-dir")
    display.DisplayMessage(gp_Pnt(0, 0, 1000), "z-dir")

def display_faces(display, shape):
    face_explorer = TopExp_Explorer(shape, TopAbs_FACE)
    face_count = 0
    while face_explorer.More():
        face = topods.Face(face_explorer.Current())
        ais_face = AIS_Shape(face)
        display.Context.Display(ais_face, False)
        label = f"Face {face_count}"
        
        bbox = Bnd_Box()
        brepbndlib.Add(face, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        
        center = gp_Pnt((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
        
        display.DisplayMessage(center, label)
        
        face_count += 1
        face_explorer.Next()
    return face_count

def highlight_face(display, shape, face_number, total_faces):
    display.Context.RemoveAll(False)
    display_origin(display)
    display_axis_lines(display)
    
    face_explorer = TopExp_Explorer(shape, TopAbs_FACE)
    face_count = 0
    while face_explorer.More():
        face = topods.Face(face_explorer.Current())
        ais_face = AIS_Shape(face)
        if face_count == face_number:
            ais_face.SetColor(Quantity_Color(1.0, 0.0, 0.0, Quantity_TOC_RGB))
            label = f"Face {face_count} (red)"
            print(f"Face {face_count} is red")
        else:
            label = f"Face {face_count}"
        
        display.Context.Display(ais_face, False)
        
        bbox = Bnd_Box()
        brepbndlib.Add(face, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        
        center = gp_Pnt((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
        
        display.DisplayMessage(center, label)
        
        face_count += 1
        face_explorer.Next()
    
    display.FitAll()

def main():
    file_path = 'data/Pr81502 - Hub OD240mm stp.STEP'
    
    display, start_display, add_menu, add_function_to_menu = init_display()
    
    shape = read_step_file(file_path)
    if shape is None:
        return
    
    display_origin(display)
    display_axis_lines(display)
    
    total_faces = display_faces(display, shape)
    print(f"Total number of entities: {total_faces}")
    
    current_face = 0
    while True:
        user_input = input(f"Enter the face number (current: {current_face}): ").strip()
        if user_input == "":
            current_face = (current_face + 1) % total_faces
        else:
            try:
                current_face = int(user_input)
                if current_face >= total_faces or current_face < 0:
                    print(f"Invalid face number. Please enter a number between 0 and {total_faces - 1}.")
                    continue
            except ValueError:
                print("Invalid input. Please enter a valid face number or press Enter to increment.")
                continue
        
        highlight_face(display, shape, current_face, total_faces)
    
    start_display()

if __name__ == "__main__":
    main()
