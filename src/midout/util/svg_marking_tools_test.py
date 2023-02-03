from midout.util.svg_marking_tools import compare_id, extract_points_from_path_d


def test_compare_id():
    assert compare_id({'a':1, 'b':2, 'c':3}, a=[1,2], b=2) == True
    assert compare_id({'a':1, 'b':2, 'c':3}, a=1, b=3) == False
    assert compare_id({'a':1, 'b':2, 'c':3}, a=1, b=[3,4]) == False
    assert compare_id({'a':1, 'b':2, 'c':3}, d=3) == False

    assert compare_id({'coords':(0,1,2)}, coords=None) == True
    assert compare_id({'coords':(0,1,2)}, coords=(0,1,2)) == True
    assert compare_id({'coords':(0,1,2)}, coords=([0,1], None, 2)) == True
    assert compare_id({'coords':(0,1,2)}, coords=(1, None, 2)) == False
    assert compare_id({'coords':(0,1,2)}, coords=(0, None, [3,4])) == False


def test_extract_points_from_path_d():
    example_string = "M805.952,38.6274 L851.207,38.6274 L873.834,61.2548 L851.207,83.8822 C842.156 70.3058,842.156 70.3058,828.579 61.2548 C819.528 47.6784,819.528 47.6784,805.952 38.6274"
    expected_coords = [(805.952,38.6274), (851.207,38.6274),(873.834,61.2548),(851.207,83.8822), (828.579, 61.2548), (805.952, 38.6274)]
    assert extract_points_from_path_d(example_string) == expected_coords

    example_string = "M896.461,38.6274 C894.199 27.3137, 885.148 18.2627, 873.834 16 C876.097 27.3137, 885.148 36.3647, 896.461 38.6274"
    expected_coords = [(896.461,38.6274), (873.834, 16.0), (896.461, 38.6274)]
    assert extract_points_from_path_d(example_string) == expected_coords

