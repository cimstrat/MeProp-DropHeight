from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopoDS import topods
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Display.SimpleGui import init_display
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.AIS import AIS_Shape
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.gp import gp_Pnt, gp_Ax1, gp_Trsf, gp_Dir
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeEdge
from OCC.Core.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.Geom import Geom_Line
from OCC.Core.BRep import BRep_Tool
from math import cos, sin, radians, pi, sqrt

from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.Geom import Geom_BSplineSurface, Geom_Plane, Geom_CylindricalSurface, Geom_ConicalSurface, Geom_SphericalSurface
from OCC.Core.GeomConvert import geomconvert_SurfaceToBSplineSurface


def read_step_file(file_path):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(file_path)
    if status != 1:
        print("Error: Cannot read STEP file")
        return None
    step_reader.TransferRoots()
    shape = step_reader.Shape()
    return shape

def read_data_file(data_file):
    blade_face, hub_top, hub_bottom, height, theta, nbr_circles, radii = None, None, None, None, None, None, []
    nbr_spurs, phi_angles = None, []
    
    with open(data_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            if blade_face is None:
                blade_face = int(line.strip())
            elif hub_top is None:
                hub_top = int(line.strip())
            elif hub_bottom is None:
                hub_bottom = int(line.strip())
            elif height is None:
                height = float(line.strip())
            elif theta is None:
                theta = float(line.strip())
            elif nbr_circles is None:
                nbr_circles_data = line.strip().split(',')
                nbr_circles = int(nbr_circles_data[0].strip())
                radii = [float(r.strip()) for r in nbr_circles_data[1:]]
            elif nbr_spurs is None:
                spurs_data = line.strip().split(',')
                nbr_spurs = int(spurs_data[0].strip())
                phi_angles = [float(phi.strip()) for phi in spurs_data[1:]]
    
    return blade_face, hub_top, hub_bottom, height, theta, nbr_circles, radii, nbr_spurs, phi_angles

def display_origin(display):
    origin = gp_Pnt(0, 0, 0)
    vertex_O = BRepBuilderAPI_MakeVertex(origin).Vertex()
    ais_vertex_O = AIS_Shape(vertex_O)
    display.Context.Display(ais_vertex_O, False)
    display.DisplayMessage(origin, "O")
    return origin

def highlight_faces(display, shape, blade_face_num, hub_top_num, hub_bottom_num):
    display.Context.RemoveAll(False)
    origin = display_origin(display)
    
    face_explorer = TopExp_Explorer(shape, TopAbs_FACE)
    face_count = 0
    blade_face, hub_top_face, hub_bottom_face = None, None, None  # Initialize variables
    while face_explorer.More():
        face = topods.Face(face_explorer.Current())
        ais_face = AIS_Shape(face)
        label = None  # Initialize label
        
        if face_count == blade_face_num:
            ais_face.SetColor(Quantity_Color(1.0, 0.0, 0.0, Quantity_TOC_RGB))  # Red
            label = "blade_face"
            blade_face = face  # Capture the blade face
            display.Context.Display(ais_face, False)
        elif face_count == hub_top_num:
            ais_face.SetColor(Quantity_Color(0.0, 1.0, 0.0, Quantity_TOC_RGB))  # Green
            label = "hub_top"
            hub_top_face = face  # Capture the hub_top face
            display.Context.Display(ais_face, False)
        elif face_count == hub_bottom_num:
            ais_face.SetColor(Quantity_Color(0.0, 0.0, 1.0, Quantity_TOC_RGB))  # ?
            label = "hub_bottom"
            hub_bottom_face = face  # Capture the hub_top face
            display.Context.Display(ais_face, False)
        else:
            display.Context.Erase(ais_face, False)  # Hide other faces

        if label is not None:
            bbox = Bnd_Box()
            brepbndlib.Add(face, bbox)
            xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
            
            center = gp_Pnt((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
            
            display.DisplayMessage(center, label)
        
        face_count += 1
        face_explorer.Next()
    
    display.FitAll()
    return origin, blade_face, hub_top_face, hub_bottom_face  # Return three faces

def calculate_and_print_distance_from_origin_to_hub_top(origin, hub_top_face):
    dist_shape_shape = BRepExtrema_DistShapeShape()
    dist_shape_shape.LoadS1(hub_top_face)
    dist_shape_shape.LoadS2(BRepBuilderAPI_MakeVertex(origin).Shape())
    dist_shape_shape.Perform()
    if dist_shape_shape.IsDone() and dist_shape_shape.Value() < float('inf'):
        distance = dist_shape_shape.Value()
        print(f"Distance from origin to hub_top surface: {distance:.3f} units")
    else:
        print("Failed to calculate distance or no valid distance found.")

def calculate_and_print_distance_from_origin_to_hub_bottom(origin, hub_bottom_face):
    dist_shape_shape = BRepExtrema_DistShapeShape()
    dist_shape_shape.LoadS1(hub_bottom_face)
    dist_shape_shape.LoadS2(BRepBuilderAPI_MakeVertex(origin).Shape())
    dist_shape_shape.Perform()
    if dist_shape_shape.IsDone() and dist_shape_shape.Value() < float('inf'):
        distance = dist_shape_shape.Value()
        print(f"Distance from origin to hub_bottom surface: {distance:.3f} units")
    else:
        print("Failed to calculate distance or no valid distance found.")

def create_point_c(display, origin, height):
    point_c = gp_Pnt(height, 0, 0)
    vertex_C = BRepBuilderAPI_MakeVertex(point_c).Vertex()
    ais_vertex_C = AIS_Shape(vertex_C)
    display.Context.Display(ais_vertex_C, False)
    display.DisplayMessage(point_c, "C")
    
    edge_OC = BRepBuilderAPI_MakeEdge(origin, point_c).Edge()
    ais_edge_OC = AIS_Shape(edge_OC)
    ais_edge_OC.SetColor(Quantity_Color(0.0, 0.0, 1.0, Quantity_TOC_RGB))  # Blue
    display.Context.Display(ais_edge_OC, False)
    
    display.FitAll()
    return point_c

def create_point_a1(display, origin):
    point_a1 = gp_Pnt(210, 0, 0)
    vertex_A1 = BRepBuilderAPI_MakeVertex(point_a1).Vertex()
    ais_vertex_A1 = AIS_Shape(vertex_A1)
    display.Context.Display(ais_vertex_A1, False)
    display.DisplayMessage(point_a1, "A1")
    
    display.FitAll()

def create_point_a2(display, origin):
    point_a2 = gp_Pnt(-135, 0, 0)
    vertex_A2 = BRepBuilderAPI_MakeVertex(point_a2).Vertex()
    ais_vertex_A2 = AIS_Shape(vertex_A2)
    display.Context.Display(ais_vertex_A2, False)
    display.DisplayMessage(point_a2, "A2")
    
    display.FitAll()

def create_centre_line(display, point_c, blade_face, theta):
    # Calculate the length of the blade part
    bbox = Bnd_Box()
    brepbndlib.Add(blade_face, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    
    length = ymax - ymin  # Length of the blade part
    
    theta_rad = radians(theta)
    
    # Define centre_line with angle theta from the y-axis in the y-z plane
    centre_line_end = gp_Pnt(point_c.X(), point_c.Y() + length * cos(theta_rad), point_c.Z() + length * sin(theta_rad))
    edge_centre_line = BRepBuilderAPI_MakeEdge(point_c, centre_line_end).Edge()
    ais_edge_centre_line = AIS_Shape(edge_centre_line)
    ais_edge_centre_line.SetColor(Quantity_Color(0.0, 0.0, 1.0, Quantity_TOC_RGB))  # Blue
    display.Context.Display(ais_edge_centre_line, False)
    display.DisplayMessage(centre_line_end, "centre_line")
    
    display.FitAll()
    return length

def create_demi_circles(display, point_c, nbr_circles, radii, theta):
    theta_rad = radians(theta)
    for i in range(nbr_circles):
        radius = radii[i]
        
        # Define the start, mid, and end points for the arc
        start_point = gp_Pnt(
            point_c.X(),
            point_c.Y() + radius * cos(-pi / 2),
            point_c.Z() + radius * sin(-pi / 2)
        )
        
        mid_point = gp_Pnt(
            point_c.X(),
            point_c.Y() + radius,
            point_c.Z()
        )
        
        end_point = gp_Pnt(
            point_c.X(),
            point_c.Y() + radius * cos(pi / 2),
            point_c.Z() + radius * sin(pi / 2)
        )
        
        # Rotate the start, mid, and end points around the center point C
        transform = gp_Trsf()
        transform.SetRotation(gp_Ax1(point_c, gp_Dir(1, 0, 0)), theta_rad)
        
        start_point.Transform(transform)
        mid_point.Transform(transform)
        end_point.Transform(transform)
        
        # Create the arc using the transformed points
        arc = GC_MakeArcOfCircle(start_point, mid_point, end_point).Value()
        edge_arc = BRepBuilderAPI_MakeEdge(arc).Edge()
        
        ais_arc = AIS_Shape(edge_arc)
        ais_arc.SetColor(Quantity_Color(0.0, 0.0, 1.0, Quantity_TOC_RGB))  # Blue
        display.Context.Display(ais_arc, False)
        display.DisplayMessage(mid_point, f"C{i+1}")
    
    display.FitAll()

def create_spurs(display, point_c, nbr_spurs, phi_angles, length, theta):
    theta_rad = radians(theta)
    for i in range(nbr_spurs):
        phi_rad = radians(phi_angles[i])
        
        spur_end = gp_Pnt(
            point_c.X(),
            point_c.Y() + length * cos(phi_rad),
            point_c.Z() + length * sin(phi_rad)
        )
        
        # Rotate the end point around the center point C
        transform = gp_Trsf()
        transform.SetRotation(gp_Ax1(point_c, gp_Dir(1, 0, 0)), theta_rad)
        
        spur_end.Transform(transform)
        
        edge_spur = BRepBuilderAPI_MakeEdge(point_c, spur_end).Edge()
        ais_spur = AIS_Shape(edge_spur)
        ais_spur.SetColor(Quantity_Color(0.5, 0.0, 0.5, Quantity_TOC_RGB))  # Purple
        display.Context.Display(ais_spur, False)
        display.DisplayMessage(spur_end, f"Spur{i+1}")
    
    display.FitAll()

def compute_intersection_points(point_c, nbr_circles, radii, nbr_spurs, phi_angles, length, theta):
    theta_rad = radians(theta)
    points = []

    for i in range(nbr_circles):
        radius = radii[i]

        for j in range(nbr_spurs):
            phi_rad = radians(phi_angles[j])
            
            # Calculate spur end point
            spur_end = gp_Pnt(
                point_c.X(),
                point_c.Y() + length * cos(phi_rad),
                point_c.Z() + length * sin(phi_rad)
            )
            
            transform = gp_Trsf()
            transform.SetRotation(gp_Ax1(point_c, gp_Dir(1, 0, 0)), theta_rad)
            spur_end.Transform(transform)
            
            # Calculate intersection points
            dx, dy = spur_end.Y() - point_c.Y(), spur_end.Z() - point_c.Z()
            A = dx * dx + dy * dy
            B = 2 * (dx * point_c.Y() + dy * point_c.Z())
            C = point_c.Y() * point_c.Y() + point_c.Z() * point_c.Z() - radius * radius
            det = B * B - 4 * A * C
            
            if det >= 0:
                t1 = (-B + sqrt(det)) / (2 * A)
                t2 = (-B - sqrt(det)) / (2 * A)
                for t in [t1, t2]:
                    if 0 <= t <= 1:
                        x = point_c.X()
                        y = point_c.Y() + t * dx
                        z = point_c.Z() + t * dy
                        points.append((i + 1, radius, j + 1, phi_angles[j], x, y, z))
    
    return points

def display_points_and_lines(display, points, height, blade_face):
    for m, radius, n, angle, x, y, z in points:
        point = gp_Pnt(x, y, z)
        vertex = BRepBuilderAPI_MakeVertex(point).Vertex()
        ais_vertex = AIS_Shape(vertex)
        ais_vertex.SetColor(Quantity_Color(1.0, 1.0, 0.0, Quantity_TOC_RGB))  # Yellow
        display.Context.Display(ais_vertex, False)
        display.DisplayMessage(point, f"P({m},{n})")
        
        # Create line LPQ(m,n)
        end_point = gp_Pnt(x - 2 * height, y, z)
        edge_LPQ = BRepBuilderAPI_MakeEdge(point, end_point).Edge()
        ais_edge_LPQ = AIS_Shape(edge_LPQ)
        ais_edge_LPQ.SetColor(Quantity_Color(0.0, 1.0, 1.0, Quantity_TOC_RGB))  # Cyan
        display.Context.Display(ais_edge_LPQ, False)
        display.DisplayMessage(end_point, f"LPQ({m},{n})")
        
        # Create Geom_Line for intersection calculation
        geom_line = Geom_Line(point, gp_Dir(-1, 0, 0))
        
        # Find intersections
        intersections = find_intersections(geom_line, blade_face)
        
        for i, intersection in enumerate(intersections):
            vertex_Q = BRepBuilderAPI_MakeVertex(intersection).Vertex()
            ais_vertex_Q = AIS_Shape(vertex_Q)
            ais_vertex_Q.SetColor(Quantity_Color(0.0, 1.0, 0.0, Quantity_TOC_RGB))  # Green
            display.Context.Display(ais_vertex_Q, False)
            display.DisplayMessage(intersection, f"Q({m},{n}){i+1}")
    
    display.FitAll()

def find_intersections(line, surface):
    intersections = []
    geom_surface = BRep_Tool.Surface(surface)

    bspline_surface = geomconvert_SurfaceToBSplineSurface(geom_surface)
    if bspline_surface:
        print(f"Face converted to B-spline/NURBS surface")
        print(bspline_surface)
        intersector = GeomAPI_IntCS(line, bspline_surface)
    else:
        print(f"Face is another type of surface: {type(surface).__name__}")
        intersector = GeomAPI_IntCS(line, geom_surface)

#    print (geom_surface)
#    if geom_surface:
#        surf = geom_surface
#        if isinstance(surf, Geom_BSplineSurface):
#            print(f'Surface type is B-spline')
#        elif isinstance(surf, Geom_Plane):
#            print(f'Surface type is plane')
#        elif isinstance(surf, Geom_CylindricalSurface):
#            print(f'Surface_type is Cylindrical')
#        elif isinstance(surf, Geom_ConicalSurface):
#            print(f'Surface_type is Conical')
#        elif isinstance(surf, Geom_SphericalSurface):
#            print(f'Surface_type is Spherical')
#        else:
#            print(f'Surface type is Other')
#    else:
#        print(f'Surface type is unknown.')

#    intersector = GeomAPI_IntCS(line, geom_surface)
    
    if intersector.IsDone() and intersector.NbPoints() > 0:
        for i in range(intersector.NbPoints()):
            intersections.append(intersector.Point(i + 1))
    
    return intersections

def compute_distances(points, blade_face):
    result = []

    for m, radius, n, angle, x, y, z in points:
        point_P = gp_Pnt(x, y, z)
        geom_line = Geom_Line(point_P, gp_Dir(-1, 0, 0))
        intersections = find_intersections(geom_line, blade_face)
        
        distances = []
        for intersection in intersections:
            distance = point_P.Distance(intersection)
            distances.append((intersection, distance))
        result.append((m, radius, n, angle, point_P, distances))
    
    return result

def print_table(distances):
    header = f"{'arc number':<12}{'radius (mm)':<12}{'spur number':<12}{'angle (deg)':<12}{'P(m,n)':<25}{'Q(m,n)1':<25}{'dist(m,n)1':<15}{'Q(m,n)2':<25}{'dist(m,n)2':<15}"
    table = [header]
    separator = '-' * len(header)
    table.append(separator)
    for m, radius, n, angle, point_P, dists in distances:
        p_coords = f"({point_P.X():.2f}, {point_P.Y():.2f}, {point_P.Z():.2f})"
        q1_coords = f"({dists[0][0].X():.2f}, {dists[0][0].Y():.2f}, {dists[0][0].Z():.2f})" if len(dists) > 0 else ""
        dist1 = f"{dists[0][1]:.2f}" if len(dists) > 0 else ""
        q2_coords = f"({dists[1][0].X():.2f}, {dists[1][0].Y():.2f}, {dists[1][0].Z():.2f})" if len(dists) > 1 else ""
        dist2 = f"{dists[1][1]:.2f}" if len(dists) > 1 else ""
        table.append(f"{m:<12}{radius:<12}{n:<12}{angle:<12}{p_coords:<25}{q1_coords:<25}{dist1:<15}{q2_coords:<25}{dist2:<15}")
    
    for line in table:
        print(line)
    
    with open('DroProp_results.txt', 'w') as file:
        for line in table:
            file.write(line + '\n')

def main():

    file_path = 'data/Pr81502 - Hub OD240mm stp.STEP'
    data_file = 'DroProp_data.txt'
    
    display, start_display, add_menu, add_function_to_menu = init_display()
    
    shape = read_step_file(file_path)
    if shape is None:
        return
    
    blade_face_num, hub_top_num, hub_bottom_num, height, theta, nbr_circles, radii, nbr_spurs, phi_angles = read_data_file(data_file)
    
    origin, blade_face, hub_top_face, hub_bottom_face = highlight_faces(display, shape, blade_face_num, hub_top_num, hub_bottom_num)  # Updated to receive hub_top_face
    
    calculate_and_print_distance_from_origin_to_hub_top(origin, hub_top_face)  # Calculate and print distance
    calculate_and_print_distance_from_origin_to_hub_bottom(origin, hub_bottom_face)  # Calculate and print distance

   
#    origin = display_origin(display)
#    origin, blade_face = highlight_faces(display, shape, blade_face_num, hub_top_num)

    point_c = create_point_c(display, origin, height)
    create_point_a1(display, origin)
    create_point_a2(display, origin)
#    create_point_c2(display, origin)
    length = create_centre_line(display, point_c, blade_face, theta)
    create_demi_circles(display, point_c, nbr_circles, radii, theta)
    create_spurs(display, point_c, nbr_spurs, phi_angles, length, theta)
    
    intersection_points = compute_intersection_points(point_c, nbr_circles, radii, nbr_spurs, phi_angles, length, theta)
    display_points_and_lines(display, intersection_points, height, blade_face)
    
    distances = compute_distances(intersection_points, blade_face)
    print_table(distances)
    
    display.View_Iso()  # Set the initial view to isometric
    display.FitAll()
    
    # Set the view to the y-z plane normal to the x-axis
    display.View.SetProj(1, 0, 0)  # View along the x-axis
    
    start_display()

if __name__ == "__main__":
    main()
