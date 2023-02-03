import re

from  xml.etree import ElementTree

import numpy as np

SVG_NAMESPACE = "{http://www.w3.org/2000/svg}"

def parse_id(e):
    """parse ids from stim SVG elements into a helpful dictionary"""
    split = e.attrib['id'].split(':')
    out = {
        'tag': e.tag.split('}')[-1]
    }
    if split[0] == "slice":
        out['name'] = 'slice'
        out['detector'] = split[1]
        out['coords'] = tuple(float(x) for x in split[2].split('_')) if len(split)==4 else None
        out['tick'] = int(split[-1])
    elif split[0] == "qubit_dots":
        out['name'] = 'qubit_dots'
    elif split[0] == "tick_borders":
        out['name'] = 'tick_borders'
    else:
        raise ValueError(f"Not Recognised: {e.attrib['id']}")
    return out


def marker_id_string(id_dict):
    """assemble an id string for a new marker from a slice id dictionary"""
    if 'coords' in id_dict:
        cs = f":{'_'.join([str(f) for f in id_dict['coords']])}"
    else:
        cs = ''
    return f"marker:{id_dict['detector']}{cs}:{id_dict['tick']}"


def compare_id(id, coords=None, **kwargs):
    """check each kwarg is in the dictionary with the right value

    single value kwargs must match exactly
    lists are intepreted as needing to match any one element in the list

    coords is special cased, because this is a tuple or list of some length
    """

    for k,v in kwargs.items():
        if k not in id:
            return False
        if isinstance(v, list):
            if id[k] not in v:
                return False
        else:
            if id[k] != v:
                return False

    if coords is not None:
        if 'coords' not in id:
            return False
        if len(coords) != len(id['coords']):
            raise ValueError(f"coords filter length didn't match coords length: {coords} != {id['coords']}")
        for i, d in enumerate(coords):
            if d is None:
                continue
            elif isinstance(d, list):
                if id['coords'][i] not in d:
                    return False
            else:
                if id['coords'][i] != d:
                    return False

    return True

def extract_points_from_path_d(path_d_string):
    out = []
    for s in re.split(" (?=[a-zA-Z])", path_d_string):
        if s[0] in ["M", "L"]:
            x,y = s[1:].split(',')
            out.append((float(x), float(y)))
        elif s[0] == "C":
            points = re.split(", ?", s[1:])
            # you can chain these curves: every 3rd point is a real point
            for p in points[2::3]:
                x,y = p.split(' ')
                out.append((float(x), float(y)))
        else:
            raise ValueError(f"Unexpected path formatting: {s} in {path_d_string}")
    return out

def shoelace_area(points):
    """return the signed shoelace area

    assumes the points are ordered around the perimeter of the shape,
    and that the polygon is not self-overlapping
    """
    return 0.5 * np.sum(
        [
            points[i - 1][0] * points[i][1] - points[i][0] * points[i - 1][1]
            for i in range(len(points))
        ]
    )

def centroid(points):
    """returns the centroid of a 2D polygon.

    assumes the points are ordered around the perimeter of the shape,
    and that the polygon is not self-overlapping
    """
    # so, this is a bit batshit, but:
    # https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    if len(points) == 1:
        return points[0]
    A = shoelace_area(points)
    if np.isclose(A, 0):
        # the area is 0 if there are 2 points, or the points are colinear
        return np.mean(points, axis=0)
    lx = [
        (points[i - 1][0] + points[i][0])
        * (points[i - 1][0] * points[i][1] - points[i][0] * points[i - 1][1])
        for i in range(len(points))
    ]
    cx = 1 / 6.0 / A * np.sum(lx)
    ly = [
        (points[i - 1][1] + points[i][1])
        * (points[i - 1][0] * points[i][1] - points[i][0] * points[i - 1][1])
        for i in range(len(points))
    ]
    cy = 1 / 6.0 / A * np.sum(ly)
    return [cx, cy]

def marker_factory(slice_group_element, id_str, marker='o'):
    """ add a marker to the slice group element"""

    # get the surrounding path, compute centroid
    slice_path_element = slice_group_element[0]
    if slice_path_element.tag == SVG_NAMESPACE+'path':
        if 'd' not in slice_path_element.attrib:
            raise ValueError(f"Unexpected slice_path_element: {slice_path_element} : {slice_path_element.attrib}")
        all_corners = extract_points_from_path_d(slice_path_element.attrib['d'])
        # because we are geniuses and don't use closepath, remove the last corner, which is a duplicate
        all_corners = all_corners[:-1]
        c = centroid(all_corners)
    elif slice_path_element.tag == SVG_NAMESPACE+'circle':
        c = (float(slice_path_element.attrib['cx']), float(slice_path_element.attrib['cy']))
    else:
        raise ValueError(f"Unexpected slice_path_element.tag={slice_path_element.tag}")

    if marker == 'o':
        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'circle')
        marker.set('r', str(3.5))
        marker.set('cx', str(c[0]))
        marker.set('cy', str(c[1]))
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        if id_str:
            marker.set('id', id_str)
    elif marker == '-':
        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'path')
        marker.set('d', f"M{c[0]-3.5},{c[1]} L{c[0]+3.5},{c[1]}")
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        if id_str:
            marker.set('id', id_str)

    elif marker == '+':
        a = 3.5

        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'path')
        marker.set(
            'd',
            f"M{c[0]-a},{c[1]} "
            f"L{c[0]+a},{c[1]} "
            f"M{c[0]},{c[1]-a} "
            f"L{c[0]},{c[1]+a}"
        )
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        if id_str:
            marker.set('id', id_str)

    elif marker == 'x':
        a = 3.0

        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'path')
        marker.set(
            'd',
            f"M{c[0]-a},{c[1]-a} "
            f"L{c[0]+a},{c[1]+a} "
            f"M{c[0]+a},{c[1]-a} "
            f"L{c[0]-a},{c[1]+a}"
        )
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        marker.set('id', id_str)

    elif marker == '*':
        a = 3.5
        b = 2.5 # approx 3.5/sqrt(2)

        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'path')
        marker.set(
            'd',
            f"M{c[0]-a},{c[1]} "
            f"L{c[0]+a},{c[1]} "
            f"M{c[0]},{c[1]-a} "
            f"L{c[0]},{c[1]+a} "
            f"M{c[0]-b},{c[1]-b} "
            f"L{c[0]+b},{c[1]+b} "
            f"M{c[0]+b},{c[1]-b} "
            f"L{c[0]-b},{c[1]+b}"
        )
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        marker.set('id', id_str)

    elif marker == 'b':
        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'path')
        marker.set(
            'd', 
            f"M{c[0]-2.5},{c[1]-2.5} "
            f"L{c[0]-2.5},{c[1]+2.5} "
            f"L{c[0]+2.5},{c[1]-2.5} "
            f"L{c[0]+2.5},{c[1]+2.5} "
            "Z"
        )
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        marker.set('id', id_str)

    elif marker == 's':
        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'path')
        marker.set(
            'd', 
            f"M{c[0]-2.5},{c[1]-2.5} "
            f"L{c[0]+2.5},{c[1]-2.5} "
            f"L{c[0]+2.5},{c[1]+2.5} "
            f"L{c[0]-2.5},{c[1]+2.5} "
            "Z"
        )
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        marker.set('id', id_str)

    elif marker == 'd':
        a = 3.5

        marker = ElementTree.SubElement(slice_group_element, SVG_NAMESPACE+'path')
        marker.set(
            'd',
            f"M{c[0]-a},{c[1]} "
            f"L{c[0]},{c[1]+a} "
            f"L{c[0]+a},{c[1]} "
            f"L{c[0]},{c[1]-a} "
            "Z"
        )
        marker.set('stroke', '#000000')
        marker.set('fill', 'none')
        if id_str:
            marker.set('id', id_str)

    else:
        raise ValueError(f"Unrecognised Marker: {marker}")

def mark_slices(tree, marker='x', **kwargs):
    """adds x markers to all slices compatible with kwargs
    
    see compare_id for argument specifications
    """
    c=0
    for e in tree.findall(f"./{SVG_NAMESPACE}g"):
        idd = parse_id(e)
        if compare_id(idd, tag='g', name='slice', **kwargs):
            marker_factory(e, id_str=marker_id_string(idd), marker=marker)
            c += 1
    return c

def color_slices(tree, color='#59FF7A', **kwargs):
    """changes the color of the background fill of the given slice
    
    see compare_id for argument specifications
    """
    c=0
    for e in tree.findall(f"./{SVG_NAMESPACE}g"):
        idd = parse_id(e)
        if compare_id(idd, tag='g', name='slice', **kwargs):
            # get the surrounding path, compute centroid
            slice_path_element = e[0]
            if slice_path_element.tag == SVG_NAMESPACE+'path' or slice_path_element.tag == SVG_NAMESPACE+'circle':
                slice_path_element.set('fill', color)
            else:
                raise ValueError(f"Unexpected slice_path_element.tag={slice_path_element.tag}")
            c += 1
    return c

def parse_string_to_tree(xml_string: str) -> ElementTree:
    """parses an xml string into an ElementTree
    
    necesary because ElementTree.parse unhelpfully only takes files not strings, 
    and ElementTree.fromstring returns an Element not an ElementTree
    """
    return ElementTree.ElementTree(element=ElementTree.fromstring(xml_string))