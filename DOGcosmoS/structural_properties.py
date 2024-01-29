import numpy as np
import unyt as U
from velociraptor import load as load_catalogue

def HI_disk_size(xyz, m, frac=0.9):
    """Calculates the radius enclosing the specified fraction of the HI mass.
    
    Arguments:
             xyz: array of dimension 3xN containing gas particle positions.
             m_HI: array of dimension 1xN containing HI masses.
             frac: mass fraction to be enclosed. Default is 0.9
    """
    transposed = False
    if xyz.ndim != 2:
        raise ValueError(
            "HI_disk_size: cannot calculate disk size for input with" " ndim != 2."
        )
    elif (xyz.shape[0] == 3) and (xyz.shape[1] == 3):
        raise ValueError(
            "HI_disk_size: cannot calculate disk size for input with" " shape (3, 3)."
        )
    elif xyz.shape[1] == 3:
        xyz = xyz.T
        transposed = True
    elif (xyz.shape[0] != 3) and (xyz.shape[1] != 3):
        raise ValueError(
            "HI_disk_size: coordinate array shape "
            + str(xyz.shape)
            + " invalid (one dim must be 3)."
        )
    
    r = np.sqrt(np.sum(np.power(xyz, 2), axis=0))
    rsort = np.argsort(r)
    rsort = np.argsort(r, kind="quicksort")
    r = r[rsort].to('kpc')
    m_HI = m_HI[rsort]
    mcumul = np.cumsum(m_HI) / np.sum(m_HI)
    
    for j in range(len(mcumul)):
        if mcumul[j] >= frac:
            R_last = r[j] #.to('kpc')
            break
   
    return R_last

def HI_disk_height(xyz, m_HI, R_halfmass, frac=0.9):
    """Calculates the height of the disk vertically enclosing the specified fraction of the HI mass within one halfmass radius.
    
    Arguments:
             xyz: array of dimension 3xN containing gas particle positions.
             m_HI: array of dimension 1xN containing HI masses.
             frac: mass fraction to be enclosed. Default is 0.9
    """
    rotmat = L_align(
        xyz, sg.gas.velocities, m_HI,
        frac=0.3,
        Laxis="z",
    )
    
    xyz = np.dot(rotmat.T, xyz.T).to('kpc').T
    #z = np.dot(rotmat.T, sg.gas.coordinates.T).to('kpc')[2,:]
    clip = np.argwhere(np.sqrt(np.sum(np.power(xyz[:,:2], 2), axis=0)) > R_halfmass)
    np.delete(xyz, clip)
    zsort = np.argsort(xyz[:,2], kind="quicksort")
    xyz = xyz[zsort]
    print(xyz[0,2])
    print(xyz[-1,2])
    mcumul = np.cumsum(m_HI[zsort])/np.sum(m_HI[zsort])
    for j in range(len(mcumul)):
        if mcumul[j] > frac:
            disk_height = (xyz[j,2]/2).to('kpc')
            break
        
    return disk_height
