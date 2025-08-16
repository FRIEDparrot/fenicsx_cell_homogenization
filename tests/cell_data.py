import numpy as np

def get_cube_data():
    """
    Get points and grid data for a cube structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0]
    ], dtype=np.float64)

    grid = np.array([
        [0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7],
        [4, 5], [4, 6], [5, 7], [6, 7]
    ], dtype=np.int64)

    return points, grid, 0.3

def get_octahedron_data():
    """
    Get points and grid data for an octahedron structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [1.0, 0.5, 0.5],
        [0.5, 1.0, 0.5],
        [0.5, 0.5, 1.0]
    ])

    grid = np.array([
        [0, 1], [0, 2], [0, 3], [0, 4],
        [5, 1], [5, 2], [5, 3], [5, 4]
    ])

    return points, grid, 0.3

def get_octet_truss_data():
    """
    Get points and grid data for an octet truss structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0.5, 0.5, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.5, 0.5],
        [0.5, 1.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 1.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
    ])

    grid = np.array([
        [0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6], [0, 7], [0, 8],
        [9, 5], [9, 6], [9, 7], [9, 8], [7, 1], [7, 4], [7, 10], [7, 11],
        [6, 5], [6, 7], [8, 1], [8, 2], [8, 10], [8, 12], [8, 5], [8, 7],
        [9, 10], [9, 12], [9, 13], [9, 11], [6, 4], [6, 3], [6, 11], [6, 13],
        [5, 2], [5, 3], [5, 12], [5, 13]
    ])

    return points, grid, 0.2

def get_tetrahedron_data():
    """
    Get points and grid data for a tetrahedron structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, np.sqrt(3) / 2, 0.0],
        [0.5, np.sqrt(3) / 6, np.sqrt(6) / 3]
    ])

    grid = np.array([
        [0, 1], [1, 2], [2, 0],
        [0, 3], [1, 3], [2, 3]
    ])

    return points, grid, 0.25

def get_cross_lattice_data():
    """
    Get points and grid data for a cross lattice structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 1.0]
    ])
    grid = np.array([
        [0, 1], [0, 2], [0, 3],
        [1, 2], [1, 3],
        [2, 3]
    ])

    return points, grid, 0.15

def get_star_cube_data():
    """
    Get points and grid data for a star cube structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0.0, 0.0, 0.0],     # 1
        [0.5, 0.5, 0.5],     # 2 center
        [1.0, 1.0, 1.0],     # 3
        [1.0, 0.0, 0.0],     # 4
        [0.0, 1.0, 1.0],     # 5
        [0.0, 1.0, 0.0],     # 6
        [1.0, 0.0, 1.0],     # 7
        [1.0, 1.0, 0.0],     # 8
        [0.0, 0.0, 1.0],     # 9
    ])

    grid = np.array([
        [0, 1], [1, 2], [3, 1], [1, 4],
        [5, 1], [1, 6], [7, 1], [1, 8]
    ])

    return points, grid, 0.2

def get_enhanced_voxel_data():
    """
    Get points and grid data for an enhanced voxel structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0, 0, 0],           # 1
        [0.5, 0, 0.5],       # 2
        [1, 0, 1],           # 3
        [1, 0, 0],           # 4
        [0, 0, 1],           # 5
        [0, 0.5, 0.5],       # 6
        [0, 1, 1],           # 7
        [0, 1, 0],           # 8
        [0.5, 0.5, 0],       # 9
        [1, 1, 0],           # 10
        [0.5, 0.5, 1],       # 11
        [1, 1, 1],           # 12
        [0.5, 1, 0.5],       # 13
        [1, 0.5, 0.5],       # 14
        [0.5, 0.5, 0.5],     # 15 center
    ])

    grid = np.array([
        [0, 1], [1, 2], [3, 1], [1, 4], [0, 5], [5, 6], [7, 5], [5, 4],
        [0, 8], [8, 9], [3, 8], [8, 7], [4, 10], [10, 11], [2, 10], [10, 6],
        [7, 12], [12, 11], [9, 12], [12, 6], [3, 13], [13, 11], [9, 13], [13, 2],
        [0, 3], [0, 7], [0, 4], [4, 2], [4, 6], [7, 9], [7, 6], [6, 11],
        [3, 9], [3, 2], [2, 11], [9, 11], [0, 14], [14, 11], [3, 14], [14, 6],
        [7, 14], [14, 2], [9, 14], [14, 4]
    ])

    return points, grid, 0.15

def get_tesseract_data():
    """
    Get points and grid data for a tesseract structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0, 0, 0],       # 1
        [1, 0, 0],       # 2
        [0.25, 0.25, 0.25],
        [0.75, 0.25, 0.25],
        [0, 1, 0],       # 5
        [0.25, 0.75, 0.25],
        [0, 0, 1],
        [0.25, 0.25, 0.75],
        [0.75, 0.75, 0.25],
        [0.75, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.25, 0.75, 0.75],
        [1, 1, 0],       # 13
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1]
    ], dtype=np.float64)

    grid = np.array([
        [0, 1], [2, 3], [0, 4], [2, 5], [0, 6], [2, 7], [8, 3], [8, 5],
        [8, 9], [10, 3], [10, 7], [10, 9], [11, 5], [11, 7], [11, 9],
        [0, 2], [1, 3], [12, 8], [4, 5], [6, 7], [13, 10], [14, 9],
        [15, 11], [6, 13], [6, 15], [4, 12], [4, 15], [15, 14], [1, 12],
        [1, 13], [13, 14], [12, 14]
    ], dtype=np.int64)

    return points, grid, 0.1

def get_vintiles_data():
    """
    Get points and grid data for a vintiles structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0, 0.5, 0.25], [0, 0.25, 0], [0, 0.75, 0], [0.25, 0.5, 0.5],
        [0.5, 1, 0.25], [0.5, 0.75, 0.5], [1, 0.5, 0.25], [0.75, 0.5, 0.5],
        [0.5, 0, 0.25], [0.75, 0, 0], [0.25, 0, 0], [0.5, 0.25, 0.5],
        [0, 0.5, 0.75], [0, 0.25, 1], [0, 0.75, 1], [0.5, 1, 0.75],
        [1, 0.5, 0.75], [0.5, 0, 0.75], [0.75, 0, 1], [0.25, 0, 1],
        [0.25, 1, 0], [0.75, 1, 0], [1, 0.75, 0], [1, 0.25, 0],
        [0.25, 1, 1], [0.75, 1, 1], [1, 0.75, 1], [1, 0.25, 1]
    ], dtype=np.float64)

    grid = np.array([
        [0, 1], [0, 2], [0, 3], [4, 5], [6, 7], [8, 9], [8, 10], [8, 11],
        [12, 13], [12, 14], [12, 3], [15, 5], [16, 7], [17, 18], [17, 19],
        [17, 11], [2, 20], [21, 22], [23, 9], [3, 11], [7, 11], [3, 5],
        [7, 5], [1, 10], [14, 24], [25, 26], [27, 18], [13, 19],
        [4, 21], [4, 20], [15, 25], [15, 24], [6, 23], [6, 22],
        [16, 27], [16, 26]
    ], dtype=np.int64)

    return points, grid, 0.2

def get_face_centered_cubic_data():
    """
    Get points and grid data for a face-centered cubic (FCC) structure.
    Based on the data from face_center.txt.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0.0, 0.0, 0.0],     # 1
        [0.5, 0.0, 0.5],      # 2
        [1.0, 0.0, 1.0],      # 3
        [1.0, 0.0, 0.0],      # 4
        [0.0, 0.0, 1.0],      # 5
        [0.0, 0.5, 0.5],      # 6
        [0.0, 1.0, 1.0],      # 7
        [0.0, 1.0, 0.0],      # 8
        [0.5, 0.5, 0.0],      # 9
        [1.0, 1.0, 0.0],      # 10
        [0.5, 0.5, 1.0],      # 11
        [1.0, 1.0, 1.0],      # 12
        [0.5, 1.0, 0.5],      # 13
        [1.0, 0.5, 0.5],      # 14
    ], dtype=np.float64)

    grid = np.array([
        [0, 1],   # STRUT 1
        [1, 2],   # STRUT 2
        [3, 1],   # STRUT 3
        [1, 4],   # STRUT 4
        [0, 5],   # STRUT 5
        [5, 6],   # STRUT 6
        [7, 5],   # STRUT 7
        [5, 4],   # STRUT 8
        [0, 8],   # STRUT 9
        [8, 9],   # STRUT 10
        [3, 8],   # STRUT 11
        [8, 7],   # STRUT 12
        [4, 10],  # STRUT 13
        [10, 11], # STRUT 14
        [2, 10],  # STRUT 15
        [10, 6],  # STRUT 16
        [7, 12],  # STRUT 17
        [12, 11], # STRUT 18
        [9, 12],  # STRUT 19
        [12, 6],  # STRUT 20
        [3, 13],  # STRUT 21
        [13, 11], # STRUT 22
        [9, 13],  # STRUT 23
        [13, 2],  # STRUT 24
        [0, 3],   # STRUT 25
        [0, 7],   # STRUT 26
        [0, 4],   # STRUT 27
        [4, 2],   # STRUT 28
        [4, 6],   # STRUT 29
        [7, 9],   # STRUT 30
        [7, 6],   # STRUT 31
        [6, 11],  # STRUT 32
        [3, 9],   # STRUT 33
        [3, 2],   # STRUT 34
        [2, 11],  # STRUT 35
        [9, 11]   # STRUT 36
    ], dtype=np.int64)

    return points, grid, 0.15

def get_volume_center_data():
    """
    Get points and grid data for a star structure.
    :return: tuple of (points, grid, diameter)
    """
    points = np.array([
        [0, 0, 0],       # 1
        [1, 0, 0],       # 2
        [0, 1, 0],       # 3
        [0, 0, 1],       # 4
        [0.5, 0.5, 0.5], # 5 (center)
        [1, 1, 1],       # 6
        [0, 1, 1],       # 7
        [1, 0, 1],       # 8
        [1, 1, 0],       # 9
    ], dtype=np.float64)

    grid = np.array([
        [0, 1], [0, 2], [0, 3], [0, 4], [4, 5], [1, 4], [4, 6], [2, 4],
        [4, 7], [8, 4], [4, 3], [3, 7], [3, 6], [2, 8], [2, 6], [6, 5],
        [1, 8], [1, 7], [7, 5], [8, 5]
    ], dtype=np.int32)

    return points, grid, 0.2
