import automol

geom_str = """
N       -0.1871898780      0.8370625452     -0.0162907117                 
O       -1.4279904629      0.4270263866      0.1251040690                 
O        0.7107966699      0.0554790352      0.0541890999                 
H        0.0099231600      1.8342759450     -0.1792209860"""

geom = automol.geom.from_string(geom_str)
zmat = automol.geom.zmatrix(geom)
zmat_str = automol.zmatrix.string(zmat)
print(zmat)
