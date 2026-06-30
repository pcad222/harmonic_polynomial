import numpy as np
import pandas as pd
import os, warnings
import matplotlib.pyplot as plt
from numpy.linalg import lstsq
from pathlib import Path

L_MAX = 3

def build_pi_basis_matrix(df, l_max = L_MAX):
    """
    Constructs harmonic polynomial pi basis  for Bx, By, and Bz up to degree l_max, finds the Glm coeff
    minimizing Bx, By and Bz simultaneouly using lstsq (np.linsalg)
    
    Parameters: 
        df :pd.DataFrame of columns ['x', 'y', 'z'] in cartesian  coordinates.
        l_max : int, takes 0-3; for values > 3,  function will take default l_max=3,raising UserWarning.

    Returns:
        basis_df(symbols only) : pd.DataFrame of columns ['l', 'm', 'pi_x', 'pi_y', 'pi_z']
        glm_coeff_df: pd.DataFrame of columns ['l', 'm', G_lm]
        
    Notes: g_lm coeff with magnitude <1e-6 are set to 0         """

    if l_max > 3:
        warnings.warn(" ***This functions only supports to l=3, so takes l=3 as default    ***",UserWarning)
       
        l_max=3

    x, y, z = [df[col_name].to_numpy() for col_name in ('x', 'y', 'z')]

    def harmonic_basis_Bx(x, y, z, l_max):
        l_list, m_list, label_list, basis_list = [], [], [], []
        for l in range(l_max + 1):
            m = np.arange(-l - 1, l + 2)
            l_list += [l] * len(m)
            m_list += list(m)
            if l == 0:
                basis_list += [0*x, 0*x, 1 + 0*x]
                label_list += ['0', '0', '1']
            elif l == 1:
                basis_list += [y, 0*y, -0.5*x, z, x]
                label_list += ['y', '0', '-1/2 x', 'z', 'x']
            elif l == 2:
                basis_list += [2*x*y, 2*y*z, -0.5*x*y, -x*z,
                               -0.25*(3*x**2 + y**2 - 4*z**2), 2*x*z, x**2 - y**2]
                label_list += ['2xy', '2yz', '-1/2 xy', '-xz',
                               '-1/4(3x²+y²−4z²)', '2xz', 'x²−y²']
            else:
                basis_list += [3*x**2*y - y**3, 6*x*y*z,
                               -0.5*(3*x**2*y + y**3 - 6*y*z**2), -1.5*x*y*z,
                               (3/8)*(x**3 + x*y**2 - 4*x*z**2),
                               -0.25*(9*x**2*z + 3*y**2*z - 4*z**3),
                               -x**3 + 3*x*z**2, 3*(x**2*z - y**2*z), x**3 - 3*x*y**2]
                label_list += ['3x²y−y³', '6xyz', '-1/2(3x²y+y³−6yz²)', '-3/2 xyz',
                               '3/8(x³+xy²−4xz²)', '-1/4(9x²z+3y²z−4z³)',
                               '-x³+3xz²', '3(x²z−y²z)', 'x³−3xy²']
        return l_list, m_list, label_list, basis_list

    def harmonic_basis_By(x, y, z, l_max):
        l_list, m_list, label_list, basis_list = [], [], [], []
        for l in range(l_max + 1):
            m = np.arange(-l - 1, l + 2)
            l_list += [l] * len(m)
            m_list += list(m)
            if l == 0:
                basis_list += [1 + 0*y, 0*y, 0*y]
                label_list += ['1', '0', '0']
            elif l == 1:
                basis_list += [x, z, -0.5*y, 0*y, -y]
                label_list += ['x', 'z', '-1/2 y', '0', '-y']
            elif l == 2:
                basis_list += [x**2 - y**2, 2*x*z, -0.25*(x**2 + 3*y**2 - 4*z**2),
                               -y*z, -0.5*x*y, -2*y*z, -2*x*y]
                label_list += ['x²−y²', '2xz', '-1/4(x²+3y²−4z²)', '-yz',
                               '-1/2 xy', '-2yz', '-2xy']
            else:
                basis_list += [x**3 - 3*x*y**2, 3*(x**2*z - y**2*z),
                               -0.5*(x**3 + 3*y**2*x - 6*x*z**2),
                               -0.25*(3*x**2*z + 9*y**2*z - 4*z**3),
                               (3/8)*(x**2*y + y**3 - 4*y*z**2), -1.5*x*y*z,
                               -3*y*z**2 + y**3, -6*x*y*z, -3*x**2*y + y**3]
                label_list += ['x³−3xy²', '3(x²z−y²z)', '-1/2(x³+3y²x−6xz²)',
                               '-1/4(3x²z+9y²z−4z³)', '3/8(x²y+y³−4yz²)',
                               '-3/2 xyz', '-3yz²+y³', '-6xyz', '-3x²y+y³']
        return l_list, m_list, label_list, basis_list

    def harmonic_basis_Bz(x, y, z, l_max):
        l_list, m_list, label_list, basis_list = [], [], [], []
        for l in range(l_max + 1):
            m = np.arange(-l - 1, l + 2)
            l_list += [l] * len(m)
            m_list += list(m)
            if l == 0:
                basis_list += [0*z, 1 + 0*z, 0*z]
                label_list += ['0', '1', '0']
            elif l == 1:
                basis_list += [0*z, y, z, x, 0*z]
                label_list += ['0', 'y', 'z', 'x', '0']
            elif l == 2:
                basis_list += [0*z, 2*x*y, 2*y*z, z**2 - 0.5*(x**2 + y**2),
                               2*x*z, x**2 - y**2, 0*z]
                label_list += ['0', '2xy', '2yz', 'z²−1/2(x²+y²)', '2xz', 'x²−y²', '0']
            else:
                basis_list += [0*z, 3*x**2*y - y**3, 6*x*y*z,
                               3*y*z**2 - (3/4)*(x**2*y + y**3),
                               z**3 - (3/2)*z*(x**2 + y**2),
                               3*x*z**2 - (3/4)*(x**3 + x*y**2),
                               3*(x**2*z - y**2*z), x**3 - 3*x*y**2, 0*z]
                label_list += ['0', '3x²y−y³', '6xyz',
                               '3yz²−¾(x²y+y³)', 'z³−³⁄₂z(x²+y²)',
                               '3xz²−¾(x³+xy²)', '3(x²z−y²z)', 'x³−3xy²', '0']
        return l_list, m_list, label_list, basis_list

    l_list, m_list, label_bx, basis_bx = harmonic_basis_Bx(x, y, z, l_max)
    _, _, label_by, basis_by = harmonic_basis_By(x, y, z, l_max)
    _, _, label_bz, basis_bz = harmonic_basis_Bz(x, y, z, l_max)

    basis_func_df = pd.DataFrame({'l': l_list,'m':m_list,'pi_x':label_bx,
    'pi_y':label_by, 'pi_z':label_bz})
    
    
    basis_matrix_Bx= np.column_stack(basis_bx)
    basis_matrix_By= np.column_stack(basis_by)
    basis_matrix_Bz= np.column_stack(basis_bz)
    basis_matrix= np.vstack((basis_matrix_Bx, basis_matrix_By, basis_matrix_Bz))


    return basis_func_df, basis_matrix

# use the df to compute the coefficients
def compute_glm_coefficients(df, l_max=L_MAX):
    basis_df, basis_matrix = build_pi_basis_matrix(df, l_max=l_max)

    bx, by, bz = [df[col] for col in ("bx", "by", "bz")]
    field = np.concatenate([bx, by, bz])

    glm_coeff = lstsq(basis_matrix, field, rcond=None)[0]
    glm_coeff[np.abs(glm_coeff) < 1e-6] = 0

    glm_df = basis_df[["l", "m"]].copy()
    glm_df["g_lm"] = glm_coeff

    return glm_df

# reconstruct field to compute rresidual 