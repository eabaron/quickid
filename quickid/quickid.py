######################################################################
## Filename:      quickid.py
## Author:        Eddie Baron <baron@ou.edu>
## Created at:    Mon Apr 24 09:31:07 2023
## Modified at:   Tue Jul  4 09:41:06 2023
## Modified by:   Eddie Baron <ebaron@psi.edu>
## Description:   first go code for quickid
######################################################################
import numpy as np
import pandas as pd
import astropy.constants as C
import astropy.units as U
from  mendeleev import element

pi=np.pi

def _count_generator(reader):
  b = reader(1024 * 1024)
  while b:
    yield b
    b = reader(1024 * 1024)

def count_lines(infile):
  with open(infile, 'rb') as fp:
    c_generator = _count_generator(fp.raw.read)
    # count each \n
    count = sum(buffer.count(b'\n') for buffer in c_generator)
  return count

def getend(filename):
  phrase = "***"
  with open(filename,'r') as f:
    for (i, line) in enumerate(f):
      if phrase in line:
        return i+1
  return -1

def read_stout(file,head):
  nl = count_lines(file)
  nlast = getend(file)
  nskip = nl - nlast + 1
  df_stout = pd.read_csv(file,sep="\s+",skiprows=1,skipfooter=nskip,\
                         names=head,engine='python')
  return df_stout


def read_atll(file,head):
  df_atll = pd.read_csv(file,sep="\s+",skiprows=1, \
                        names=head,engine='python')
      
  return df_atll

def flu_mihalas(Aul,gu,gl,nu):
  flu = gu/gl*C.m_e.cgs.value*C.c.cgs.value**2*Aul/ \
    (4*pi*C.alpha.value*C.h.cgs.value*nu**2)
  return flu

def gettrans_stout(DIR,sfile):


  head_tp = ["flag","levl","levu","Aij","strans"]
  head_nrg = ["level","E","g","term"]
  df_lin = read_stout(sfile+".tp",head_tp)
  df_lvl = read_stout(sfile+".nrg",head_nrg)

  # El = np.empty(len(df_tp))
  # gl = np.empty(len(df_tp))
  # Eu = np.empty(len(df_tp))
  # gu = np.empty(len(df_tp))
  # wl = np.empty(len(df_tp))
  # nu = np.empty(len(df_tp))
  # nu = np.empty(len(df_tp))

  # for i_ in range(len(df_tp)):
  #   i, = np.where(df_nrg["level"] == df_tp["levl"][i_])
  #   j, = np.where(df_nrg["level"] == df_tp["levu"][i_])
  #   El[i_] = df_nrg["E"].iloc[i[0]]
  #   gl[i_] = df_nrg["g"].iloc[i[0]]
  #   Eu[i_] = df_nrg["E"].iloc[j[0]]
  #   gu[i_] = df_nrg["g"].iloc[j[0]]

  levels = df_lvl.level.to_list()
  indl = [ levels.index(value) for value in df_lin.levl ]
  indu = [ levels.index(value) for value in df_lin.levu ]
  El = df_lvl['E'].loc[indl].to_numpy()
  gl = df_lvl['g'].loc[indl].to_numpy()
  Eu = df_lvl['E'].loc[indu].to_numpy()
  gu = df_lvl['g'].loc[indu].to_numpy()

  wl = (1.0/(Eu-El)) * U.cm
  nu = (Eu - El)*C.c.cgs.value
  wl = wl.to("Angstrom").value
  flu = flu_mihalas(df_lin["Aij"],gu,gl,nu)

  df_lin['wl'] = wl
  df_lin['E_l'] = El
  df_lin['g_l'] = gl
  df_lin['E_u'] = Eu
  df_lin['g_u'] = gu
  with np.errstate(divide='ignore'):
    df_lin['loggf'] = np.log10(df_lin['g_l']*flu)
  # remove transitions without A values
  return df_lin[(df_lin.Aij != 0.0)].copy(),df_lvl


def gettrans_atll(DIR,ion):
  ############ 1      2      3      4       5      6    7    8
  head_lin = ["wl","wl_d","Aij","Aij_d","levl","levu","ind","Z",\
              #  9      10     11
              "ion","trans","flag"]
  ##############  1      2      3      4       5      6
  head_lvl = ["level","config","iel","term","iterm","isgnd",\
              #  7      8       9      10    11   12    13
              "ismeta","char","iq1","iq2","iq3","glo","ghi",\
              # 14     15     16    17    18   19    20   21
              "par","isun","coup","nup","E","dE","flag","litref"]

  df_lin = read_atll(DIR+"plin."+ion,head_lin)
  df_lvl = read_atll(DIR+"plvl."+ion,head_lvl)



  # El = np.empty(len(df_lin))
  # gl = np.empty(len(df_lin))
  # Eu = np.empty(len(df_lin))
  # gu = np.empty(len(df_lin))
  # wl = np.empty(len(df_lin))
  # nu = np.empty(len(df_lin))
  # nu = np.empty(len(df_lin))

  # for i_ in range(len(df_lin)):
  #   i, = np.where(df_lvl["level"] == df_lin["levl"][i_])
  #   j, = np.where(df_lvl["level"] == df_lin["levu"][i_])
  #   El[i_] = df_lvl["E"].iloc[i[0]]
  #   gl[i_] = df_lvl["ghi"].iloc[i[0]]
  #   Eu[i_] = df_lvl["E"].iloc[j[0]]
  #   gu[i_] = df_lvl["ghi"].iloc[j[0]]

  levels = df_lvl.level.to_list()
  indl = [ levels.index(value) for value in df_lin.levl ]
  indu = [ levels.index(value) for value in df_lin.levu ]
  El = df_lvl['E'].loc[indl].to_numpy()
  gl = df_lvl['ghi'].loc[indl].to_numpy()
  isgnd = df_lvl['isgnd'].loc[indl].to_list()
  ismeta = df_lvl['ismeta'].loc[indl].to_list()

  Eu = df_lvl['E'].loc[indu].to_numpy()
  gu = df_lvl['ghi'].loc[indu].to_numpy()

  wl = (1.0/(Eu-El)) * U.cm
  nu = (Eu - El)*C.c.cgs.value
  wl = wl.to("Angstrom").value
  flu = flu_mihalas(df_lin["Aij"],gu,gl,nu)

  df_lin['wl'] = wl
  df_lin['E_l'] = El
  df_lin['g_l'] = gl
  df_lin['isgnd'] = isgnd
  df_lin['ismeta'] = ismeta

  df_lin['E_u'] = Eu
  df_lin['g_u'] = gu

  with np.errstate(divide='ignore'):
    df_lin['loggf'] = np.log10(df_lin['g_l']*flu)

  # remove transitions without A values
  return df_lin[(df_lin.Aij != 0.0)].copy(),df_lvl

def sel_trans(TRANS,df_lin):
  if "A".casefold() in TRANS.casefold():
    return df_lin
  container = list("AFSPMR".casefold())
  items = list(TRANS.casefold()) 
  if not all(x in container for x in items):
    raise ValueError("Transition must be in AFSPMR")
  df_tmp = df_lin.copy()
  if "F".casefold() in TRANS.casefold():
    df_tmp = df_tmp[df_tmp.trans >= 2]
  elif "P".casefold() in TRANS.casefold():
    df_tmp = df_tmp[df_tmp.trans == 0]
  elif "S".casefold() in TRANS.casefold():
    df_tmp = df_tmp[df_tmp.trans == 1]
  if "M".casefold() in TRANS.casefold():
    df_tmp = df_tmp[df_tmp.ismeta == "T"]
  elif "R".casefold() in TRANS.casefold():
    df_tmp = df_tmp[df_tmp.isgnd == "T"]
  return df_tmp.copy()

def get_atll_lines(el,istg,wl_range,ATLL_DIR):
  elmendeleev = element(el)
  Z = elmendeleev.atomic_number
  kid = Z*100+(istg-1)
  if Z < 10:
    skid = "0"+str(kid)
  else:
    skid = str(kid)

  df_atll,df_lvl = gettrans_atll(ATLL_DIR,skid)
  lines = df_atll[(df_atll['wl'] >= wl_range[0] ) & (df_atll['wl'] <= wl_range[1])]
  return lines

def get_stout_lines(el,istg,wl_range,STOUT_DIR):
  shelp = el.lower() + "_" + str(istg)
  sfile = STOUT_DIR +  el.lower() + "/" + shelp + "/" + shelp

  df_stout,df_lvl = gettrans_stout(STOUT_DIR,sfile)
  lines = df_stout[(df_stout['wl'] >= wl_range[0] ) & (df_stout['wl'] <= wl_range[1])]
  return lines

def select_gf(lines,ax,gfmin,wlscale='angstrom',isH=False):
  hseries = ["Ly","H","Pa","Br","Pf","Hu"]
  T = [r"$\alpha$",r"$\beta$",r"$\gamma$",r"$\delta$",r"$\epsilon$"]
  def txt(levl,levu):
    if isH:
      dl = levu-levl
      if levl < len(hseries) and dl <= len(T):
        mark = hseries[levl-1]+T[dl-1]
      elif levl < len(hseries):
        mark = hseries[levl-1]+str(dl)
      else: 
        mark = str(levl) + r"$-$" + str(levu)
    else:
        mark = str(levl) + r"$-$" + str(levu)
    return mark

  if wlscale == 'angstrom':
    wfac = 1.
  elif wlscale == 'micron':
    wfac = 1e-4
  
  ymin,ymax = ax.get_ylim()
  yoffset = (ymax-ymin)/10.

  i = 0
  for line in lines.itertuples(index=False,name="Pandas"):
    if line.loggf > gfmin:
      i += 1
      ax.axvline(x=line.wl*wfac)
      ymark = ymin +  yoffset + i%2*yoffset
      ax.text(line.wl*wfac,ymark,txt(line.levl,line.levu),fontsize=6)


def from_roman(num: str) -> int:
    roman_numerals = {'I':1, 'V':5, 'X':10, 'L':50, 'C':100, 'D':500, 'M':1000}
    result = 0
    for i,c in enumerate(num):
        if (i+1) == len(num) or roman_numerals[c] >= roman_numerals[num[i+1]]:
            result += roman_numerals[c]
        else:
            result -= roman_numerals[c]
    return result

