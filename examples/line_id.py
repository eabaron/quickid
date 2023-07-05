import pandas as pd
import numpy as np
import click
from astropy import units as U
from astropy import constants as C
import quickid.quickid as QID
import pylab

pi=np.pi


def read_melissa_csv(infile,z=None):
  df = pd.read_csv(infile)

  return df

def read_smooth(infile):
  ws,fs,es = np.loadtxt(infile,unpack=True)
  return ws,fs,es


if __name__ == "__main__":
  import os.path as path
  HOME = path.expanduser("~")

  # dfile = "ACKO_MIRI_from_Chris-2023-06-27_5-13_spexsmoothed.dat"
  dfile = "ACKO_MIRI_from_Chris-13-25-2023-06-27_spexsmoothed.dat"
  fn = click.prompt("Spectrum file",default=dfile)

  DATADIR= HOME+"/projects/sn2022acko/"
  fn = DATADIR+fn
  z = click.prompt("z",default=0.0053,type=float)

  smth  = np.loadtxt(fn)

  # wl_range = [smth[0,0], smth[-1,0]]
  wl_range = [(smth[0,0] * U.micron).to("Angstrom").value, (smth[-1,0]* U.micron).to("Angstrom").value]
  # wl_range = [smth[0,0]/(1+z), smth[-1,0]/(1+z)]

  fig,ax = pylab.subplots()
  while True:
    ion = click.prompt("Give ion q/Q to quit",default="H I")
    if ion.casefold() in "Q".casefold(): break
    el = ion.split()[0].capitalize()
    istg = QID.from_roman(ion.split()[1])
    ATLL_DIR = HOME + "/projects/atll/"
    STOUT_DIR = HOME + "/projects/stout/"


    isH = False
    if el == "H":
      isH = True
    if isH:
      lines = QID.get_stout_lines(el,istg,wl_range,STOUT_DIR)
    else:
      lines = QID.get_atll_lines(el,istg,wl_range,ATLL_DIR)
  
    while True:
      # selection criteria
      if "trans" in lines:
        TRANS = click.prompt("Give TRANS (Permitted, Forbidden, Semi-Forbidden, Meta-Stable, Resonance, All",default="A")
        lines_ = QID.sel_trans(TRANS,lines)
      else:
        lines_  = lines


      GFMIN = click.prompt("Give Min log(gf) value",default=-100,type=float)
      ax.clear()
      if GFMIN >= 100: break



      p1, = ax.plot(smth[:,0],smth[:,1]/np.max(smth[:,1]),lw=1,color='black')
      # p1, = ax.plot(smth[:,0]*(1+z),smth[:,1],lw=1,color='black')

      QID.select_gf(lines_,ax,GFMIN,wlscale='micron',isH=isH)

      fig.show()

