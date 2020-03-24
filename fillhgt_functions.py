def fillhgt_falcon50000(volume, diameter=27, cone_height=14.5, cone_tipdiameter=5):
    """
    Describes 50000 ul falcons geometrically and 
    calculates the respective fill-height to a fill-volume.
    """
    # measurements: 30 x 150 mm

    cone_slope = (diameter / 2 - cone_tipdiameter / 2) / cone_height
    cone_maxvol = 0.33 * np.pi * cone_height * ((diameter / 2) ** 2 + (diameter / 2) * (cone_tipdiameter / 2) + (cone_tipdiameter / 2) ** 2)
    cone_heighttip = cone_tipdiameter / (2 * cone_slope)
    cone_voltip = 0.33 * np.pi * cone_heighttip * (cone_tipdiameter / 2) ** 2

    ##catch errors##
    if volume > 50000:
        raise ValueError("volume exceeds 50 ml falcon volume")
    if volume < 0:
        raise ValueError("negative volumes not possible")

    ##distribute volume to modeled parts##
    #volume in cone
    if volume <= cone_maxvol:
        cone_vol = volume
    else:
        cone_vol = cone_maxvol
    
    #volumein cylinder
    cyl_vol = volume - cone_vol

    ##calculate height in cone##
    # volume: 1/3 * pi * h * r^2, r = m * h
    if cone_vol == cone_maxvol:
        fillhgt_cone = cone_height # if the cone is full, then the fill-height in the cone is the cone height
    else:
        cone_fullvol = volume + cone_voltip 
        fillhgt_conefull = ((3 * cone_fullvol) / (np.pi * cone_slope ** 2)) ** 0.33
        fillhgt_cone = fillhgt_conefull - cone_heighttip

    ##top cylinder##
    if cyl_vol == 0:
        fillhgt_cyl = 0
    else:
        fillhgt_cyl = cyl_vol / (np.pi * (diameter / 2) ** 2)

    ##add heights##
    fillhgt = fillhgt_cone + fillhgt_cyl
    
    return fillhgt

#==============================================================================

def fillhgt_falcon15000(volume, diameter=14, cone_height=21, cone_tipdiameter=4):
    """
    Describes 15000 ul falcons geometrically and 
    calculates the respective fill-height to a fill-volume.
    """
    # measurements: 20 x 120 mm

    cone_slope = (diameter / 2 - cone_tipdiameter / 2) / cone_height
    cone_maxvol = 0.33 * np.pi * cone_height * ((diameter / 2) ** 2 + (diameter / 2) * (cone_tipdiameter / 2) + (cone_tipdiameter / 2) ** 2)
    cone_heighttip = cone_tipdiameter / (2 * cone_slope)
    cone_voltip = 0.33 * np.pi * cone_heighttip * (cone_tipdiameter / 2) ** 2

    ##catch errors##
    if volume > 15000:
        raise ValueError("volume exceeds 15 ml falcon volume")
    if volume < 0:
        raise ValueError("negative volumes not possible")

    ##distribute volume to modeled parts##
    #volume in cone
    if volume <= cone_maxvol:
        cone_vol = volume
    else:
        cone_vol = cone_maxvol
    
    #volumein cylinder
    cyl_vol = volume - cone_vol

    ##calculate height in cone##
    # volume: 1/3 * pi * h * r^2, r = m * h
    if cone_vol == cone_maxvol:
        fillhgt_cone = cone_height # if the cone is full, then the fill-height in the cone is the cone height
    else:
        cone_fullvol = volume + cone_voltip 
        fillhgt_conefull = ((3 * cone_fullvol) / (np.pi * cone_slope ** 2)) ** 0.33
        fillhgt_cone = fillhgt_conefull - cone_heighttip

    ##top cylinder##
    if cyl_vol == 0:
        fillhgt_cyl = 0
    else:
        fillhgt_cyl = cyl_vol / (np.pi * (diameter / 2) ** 2)

    ##add heights##
    fillhgt = fillhgt_cone + fillhgt_cyl
    
    return fillhgt