import numpy as np
import xarray as xr
from scipy.io import FortranFile

def read_var(fn, var, template, nz, nt, unit=""):
    """
    Quick function to read in Fortran arrays from the output of the
    Fortran code for the original CLUBB PDF closure.
    """
    ret = template.copy()
    variables = [
        # input
        "rtm",
        "thlm",
        "wp2",
        "wp3",
        "wprtp",
        "wpthlp",
        "rtp2",
        "thlp2",
        "rtpthlp",
        "p",
        # output
        "wp4",
        "wp2rtp",
        "wprtp2",
        "wp2thlp",
        "wpthlp2",
        "wprtpthlp",
        "wp2rcp",
        "wp2thvp",
        "wprcp",
        "wpthvp",
        "rtprcp",
        "rtpthvp",
        "thlprcp",
        "thlpthvp",
        "rcp2",
        "rtp2_clipped",
        "thlp2_clipped",
        "rtpthlp_clipped",
        "uf",
        "w_1",
        "w_2",
        "var_w_1",
        "var_w_2",
        "w_1_n",
        "w_2_n",
        "sig_w_sqd",
        "rt_1",
        "rt_2",
        "var_rt_1",
        "var_rt_2",
        "thl_1",
        "thl_2",
        "var_thl_1",
        "var_thl_2",
        "rrtthl",
        "rcm",
        "rc_1",
        "rc_2",
        "rcmi",
        "rc_1i",
        "rc_2i",
        "rsatl_1",
        "rsatl_2",
        "cf",
        "cf1",
        "cf2",
        "cfi",
        "cf1i",
        "cf2i",
        "alpha_thl",
        "alpha_rt",
        "crt_1",
        "crt_2",
        "cthl_1",
        "cthl_2",
        "chi_1",
        "chi_2",
        "stdev_chi_1",
        "stdev_chi_2",
        "stdev_eta_1",
        "stdev_eta_2",
        "covar_chi_eta_1",
        "covar_chi_eta_2",
        "corr_chi_eta_1",
        "corr_chi_eta_2",
    ]
    if var in variables:
        f = FortranFile(fn, "r")
        for i, v in enumerate(variables):
            if v == var:
                ret.values = f.read_reals(float).reshape((nz, nt), order="F")
                print("Variable {:s} ({:2d}) Found!".format(v, i))
                break
            else:
                _ = f.read_record(float)
        f.close()
        ret.name = var.upper()
        ret.attrs["long_name"] = var.upper()
        ret.attrs["unit"] = unit
        return ret
    else:
        print("Wrong Variable Name!")
        return None

def read_var_skx(fn, var, template, nz, nt, unit=""):
    """
    Quick function to read in Fortran arrays from the output of the
    Fortran code for the modified CLUBB PDF closure with skx input.
    """
    ret = template.copy()
    variables = [
        # input
        "rtm",
        "thlm",
        "wp2",
        "wp3",
        "wprtp",
        "wpthlp",
        "rtp2",
        "thlp2",
        "rtpthlp",
        "rtp3",
        "thlp3",
        "p",
        # output
        "wp4",
        "wp2rtp",
        "wprtp2",
        "wp2thlp",
        "wpthlp2",
        "wprtpthlp",
        "wp2rcp",
        "wp2thvp",
        "wprcp",
        "wpthvp",
        "rtprcp",
        "rtpthvp",
        "thlprcp",
        "thlpthvp",
        "rcp2",
        "rtp2_clipped",
        "rtp3_clipped",
        "thlp2_clipped",
        "thlp3_clipped",
        "rtpthlp_clipped",
        "uf",
        "w_1",
        "w_2",
        "var_w_1",
        "var_w_2",
        "w_1_n",
        "w_2_n",
        "sig_w_sqd",
        "rt_1",
        "rt_2",
        "var_rt_1",
        "var_rt_2",
        "thl_1",
        "thl_2",
        "var_thl_1",
        "var_thl_2",
        "rrtthl",
        "rcm",
        "rc_1",
        "rc_2",
        "rcmi",
        "rc_1i",
        "rc_2i",
        "rsatl_1",
        "rsatl_2",
        "cf",
        "cf1",
        "cf2",
        "cfi",
        "cf1i",
        "cf2i",
        "alpha_thl",
        "alpha_rt",
        "crt_1",
        "crt_2",
        "cthl_1",
        "cthl_2",
        "chi_1",
        "chi_2",
        "stdev_chi_1",
        "stdev_chi_2",
        "stdev_eta_1",
        "stdev_eta_2",
        "covar_chi_eta_1",
        "covar_chi_eta_2",
        "corr_chi_eta_1",
        "corr_chi_eta_2",
    ]
    if var in variables:
        f = FortranFile(fn, "r")
        for i, v in enumerate(variables):
            if v == var:
                ret.values = f.read_reals(float).reshape((nz, nt), order="F")
                print("Variable {:s} ({:2d}) Found!".format(v, i))
                break
            else:
                _ = f.read_record(float)
        f.close()
        ret.name = var.upper()
        ret.attrs["long_name"] = var.upper()
        ret.attrs["unit"] = unit
        return ret
    else:
        print("Wrong Variable Name!")
        return None
