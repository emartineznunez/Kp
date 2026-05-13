# =====    Import libraries of interest    ===== #
# ---------------------------------------------- #
import os                                        #
import time                                      #
# ---------------------------------------------- #
import numpy                as np                #
# ---------------------------------------------- #
import matplotlib.pyplot    as plt               #
import matplotlib.ticker    as mticker           #
# ---------------------------------------------- #
from   scipy.optimize  import root_scalar        #
# ---------------------------------------------- #
import pubchempy       as     pcp                # to retrieve from PubChem
# ---------------------------------------------- #
import rdkit                                     # SMILES --> coordinates
from   rdkit           import Chem               #
from   rdkit.Chem      import AllChem            #
# ---------------------------------------------- #
import py3Dmol                                   # 3D visualization of molecules
# ---------------------------------------------- #
import pyscf                                     # for elec struct calculations
from   pyscf.geomopt   import geometric_solver   #
from   pyscf.hessian   import thermo             #
# ---------------------------------------------- #
import ipywidgets      as     w                  # to add buttons
# ---------------------------------------------- #
from   google.colab    import files              # to access to generated files
from   IPython.utils   import io                 # to capture output
from   IPython.display import HTML               # needed for 3D visualization
from   IPython.display import display            # needed for 3D visualization
from   IPython.display import Markdown           # needed for 3D visualization
# ---------------------------------------------- #
from   constants       import m_u,m_e,q_e,h
from   constants       import k_B,c_0,eps0,NA
from   constants       import P_o,c_o,R,hbar,a_0,Eh,Hz_au
from   constants       import NPOINTSXI,REL_XI_EQ
from   constants       import ZERO1,ZERO2,ZERO3,ZERO4
from   constants       import FONTSIZE
# ============================================== #

# ------------------------------------------------
last_fig  = None
# ------------------------------------------------

# ============================================== #
TEXT1 = '''
===================
Firstly, introduce your reaction, using the following format:

      A + 2 B -> 3 C + 2 D

Make sure to include blank spaces between the stoichiometric coefficients and the chemical species.
For example, for the reaction studied in Notebook 1, you should enter:

     N2O4 -> 2 NO2
===================

'''
# ------------------------------------------------
TEXT2 = '''
~~~~~~~~~~~~~~~~~~~
Information:

 * reactants (%i): %s
 * products  (%i): %s

 * equation for the reaction: %s
~~~~~~~~~~~~~~~~~~~

'''
# ============================================== #


# ============================================= #
#          Interaction with student(s)          #
# --------------------------------------------- #
def ask_for_float(question,ntries=3):
    count = 0
    while True:
      try:
        value = float(input(question))
        break
      except:
        count += 1
        if count == ntries:
              print("      getting data failed too many times... aborting!")
              raise Exception
        print("      something went wrong... trying again...")
    return value
# --------------------------------------------- #
def ask_for_reaction(max_nerr=3):

    nus0,molecules0 = None,None

    print(TEXT1)
    answer,nerr = "n",0
    while True:
        if nerr == max_nerr:
            print("-------------------")
            print("Too many errors... aborting...")
            print("-------------------")
            return nus0,molecules0
        try:
            reaction = input(" * insert reaction: ")
            print("")
            nus,molecules = string_to_reaction(reaction)
            if nus is None: raise Exception

            #---- Print info to make sure all is correct ----
            string = reaction_to_string(nus,molecules)
            nR = len([nu for nu in nus if nu<0])
            nP = len([nu for nu in nus if nu>0])
            sR = ", ".join([molecule for nu,molecule in zip(nus,molecules) if nu < 0])
            sP = ", ".join([molecule for nu,molecule in zip(nus,molecules) if nu > 0])
            print(TEXT2%(nR,sR,nP,sP,string))
            answer = input("Is the reaction correct? (Type 'yes' or 'no'): ")
            answer = answer.strip().lower()[0]
            if answer == "y":
               print("\nGreat! Now, you are ready to study your equilibrium!\n\n")
               return nus,molecules
            else:
               print("Ok, so let us try it again! :)\n\n")
               nerr += 1
        except:
            nerr += 1
            print("-------------------")
            print("There was some kind of problem... Let us try again!")
# ============================================= #


# ============================================= #
def reaction_to_string(nus,molecules):
    stringR = []
    stringP = []
    for nu, molecule in zip(nus,molecules):
        if nu < 0: stringR += ["%s * %s(g)"%(-nu,molecule)]
        if nu > 0: stringP += ["%s * %s(g)"%(+nu,molecule)]
    return " + ".join(stringR) + " ⇌ " + " + ".join(stringP)
# --------------------------------------------- #
def string_to_reaction(string):
    reactants, products = string.split("->")
    reactants = [i.strip() for i in reactants.split("+")]
    products  = [i.strip() for i in  products.split("+")]
    nus,molecules = [],[]
    for r in reactants:
        r = [i.strip() for i in r.split()]
        if   len(r) == 1: nu,molecule = 1,r[0]
        elif len(r) == 2: nu,molecule = r
        else            : raise Exception
        nus.append(-int(nu))
        molecules.append(molecule)
    for r in products:
        r = [i.strip() for i in r.split()]
        if   len(r) == 1: nu,molecule = 1,r[0]
        elif len(r) == 2: nu,molecule = r
        else            : raise Exception
        nus.append(+int(nu))
        molecules.append(molecule)
    return np.array(nus),np.array(molecules)
# --------------------------------------------- #
def prepare_variables(T, P, nus, n_0, xi):
    T   = np.asarray(T)
    P   = np.asarray(P)
    nus = np.asarray(nus) #, dtype=float)  # (S,)
    n_0 = np.asarray(n_0) #, dtype=float)  # (S,)
    xi  = np.asarray(xi)
    # Expand for broadcasting with xi (which can be a scalar or have shape (M, N))
    # If xi.ndim == 0 (scalar), this keeps the shape (S,) -> (S,)
    # If xi.ndim == 2, it becomes (S, 1, 1) and broadcasts to (S, M, N)
    expand = (None,)*xi.ndim
    n_0_e  = n_0[(...,) + expand]
    nus_e  = nus[(...,) + expand]
    # total sums (scalars)
    dnu    = nus.sum()
    ntot_0 = n_0.sum()
    return T, P, nus_e, n_0_e, xi, dnu, ntot_0
# ============================================= #


# ============================================= #
#        Functions for downloading FILES        #
# ============================================= #
# ---- button to download file ----
def download_file(fname):
    btn = w.Button(
            description=f"Download {fname}",
            icon="download",
            button_style="primary",
            layout=w.Layout(width='250px', height='38px', margin='0 0 0 55px')
    )
    # action when clicking
    btn.on_click(lambda _: files.download(fname))
    return btn
# --------------------------------------------- #
def pyscf_download(molecule,functional,basis,DFTGRID,which_ones=[]):

    xyz_guess,xyz_opt,output_opt,output_frq = files_of_interest(molecule,functional,basis,DFTGRID)
    print(rf"     - file(s) to download:")
    if 1 in which_ones:
         print(rf"       {xyz_opt:s}")
         btn1 = download_file(xyz_opt)   ; display(btn1)
    if 2 in which_ones:
         print(rf"       {output_opt:s}")
         btn2 = download_file(output_opt); display(btn2)
    if 3 in which_ones:
         print(rf"       {output_frq:s}")
         btn3 = download_file(output_frq); display(btn3)
# ============================================= #


# ============================================= #
# Functions for PART 1: CHEMICAL THERMODYNAMICS #
# ============================================= #
# In following functions, it is assumed that    #
# all data are provided in SI units.            #
# Moreover, P is used for total pressure; lower #
# p is used for partial pressures.              #
# ============================================= #
def limits_xi(n_0,nus):
    # maximum value of xi (calculated considering consumption of reactants)
    xi_max = min([-n_0_i/nu_i for n_0_i,nu_i in zip(n_0,nus) if nu_i < 0])
    # minimum value of xi (calculated considering consumption of products)
    xi_min = max([-n_0_i/nu_i for n_0_i,nu_i in zip(n_0,nus) if nu_i > 0])
    # Ensure numerical zeros are displayed as +0.0 for clarity
    if xi_min == -0.0: xi_min = 0.0
    if xi_max == -0.0: xi_max = 0.0
    return xi_min,xi_max
# --------------------------------------------- #
def get_DHo(T,refdata):
    '''Delta{H}^o as a function of T'''
    DHo_ref,DSo_ref,DCPo_ref,TREF = refdata
    DHo_T = DHo_ref + DCPo_ref * T * (1 - TREF/T)
    return DHo_T
# --------------------------------------------- #
def get_DGo(T,refdata):
    '''Delta{G}^o as a function of T'''
    DHo_ref,DSo_ref,DCPo_ref,TREF = refdata
    DGo_T  = DHo_ref - T * DSo_ref
    DGo_T += DCPo_ref  * ( T - TREF + T*np.log(TREF / T))
    return DGo_T
# --------------------------------------------- #
def get_Gast(T,P,nus,refdata):
    DGast = get_DGo(T,refdata) + R*T*np.log(P/P_o) * nus.sum()
    return DGast
# --------------------------------------------- #
def get_DDGmix(xi, T, P, n_0, nus):
    """
    Accepts xi as either a scalar or an array (e.g., a mesh of shape (M, N)); T and p can be scalars or broadcast-compatible arrays.
    """

    T, P, nus, n_0, xi, dnu, ntot_0 = prepare_variables(T, P, nus, n_0, xi)

    # magnitudes dependent on xi
    n_xi    = n_0    + xi * nus   # (S,) o (S,M,N)
    ntot_xi = ntot_0 + xi * dnu   # scalar or (M,N)

    # Get molar fractions
    y_xi = n_xi /ntot_xi
    y_0  = n_0  /ntot_0

    # Notice that 0*ln(0) --> 0
    maskA = (n_xi > ZERO1) & (y_xi > ZERO1)
    termA = np.zeros_like(y_xi, dtype=float)
    termA[maskA] = n_xi[maskA] * np.log(y_xi[maskA])

    maskB = (n_0  > ZERO1) & (y_0  > ZERO1)
    termB = np.zeros_like(y_0, dtype=float)
    termB[maskB] = n_0[maskB]  * np.log(y_0[maskB])

    DDGmix = R*T*np.sum(termA-termB, axis=0)
    return DDGmix
# --------------------------------------------- #
def get_G_PT(xis,P,T,n_0,nus,refdata):
    '''Gibbs free energy for a given T,p; actually it is G(xi)-G(0)'''
    Gast   = get_Gast(T,P,nus,refdata)
    DDGmix = get_DDGmix(xis,T,P,n_0,nus)
    DGtot  = Gast * xis + DDGmix
    return DGtot
# --------------------------------------------- #
def get_A_VT(xis,V,T,n_0,nus,refdata):
    '''Helmholtz free energy for a given T,V'''
    ntot_0  = n_0.sum()
    ntot_xi = ntot_0 + nus.sum() * xis
    DGo_T   = get_DGo(T,refdata)
    P_xi    = ntot_xi*R*T/V
    term_lineal = DGo_T * xis
    termA   = ntot_xi * R *T * (np.log( P_xi/P_o                  )-1)
    termA  -= ntot_0  * R *T * (np.log((P_xi/P_o)*(ntot_0/ntot_xi))-1)
    DDGmix  = get_DDGmix(xis,T,P_xi,n_0,nus)
    term_nolin  = termA+DDGmix
    DAtot   = term_lineal + term_nolin
    return DAtot, term_lineal, term_nolin, P_xi
# --------------------------------------------- #
def get_Qp_PT(n_0,nus,xi,P):
    '''Quotient of reaction (pressure)'''
    # Get partial pressures for this value of xi
    n_xi = n_0 + xi*nus
    y_xi = n_xi/n_xi.sum()
    p_xi = y_xi * P
    with np.errstate(divide='ignore', invalid='ignore', over='ignore', under='ignore'):
         Qp = np.prod([(p_xi_j/P_o)**nu_j for nu_j,p_xi_j in zip(nus,p_xi)])
    return Qp
# --------------------------------------------- #
def get_Qp_VT(n_0,nus,xi,V,T):
    '''Quotient of reaction (volume)'''
    n_xi = n_0 + xi*nus
    y_xi = n_xi/(n_xi.sum())
    P    = n_xi.sum() *R*T/V
    p_xi = y_xi * P
    with np.errstate(divide='ignore', invalid='ignore', over='ignore', under='ignore'):
         Qp = np.prod([(p_xi_j/P_o)**nu_j for nu_j,p_xi_j in zip(nus,p_xi)])
    return Qp
# --------------------------------------------- #
def get_xieq_PT(P,T,n_0,nus,refdata):
    '''Extent of reaction in terms of P and T'''
    # Get Eq constant
    dGo_T = get_DGo(T,refdata)
    Kp_T  = np.exp(-dGo_T/R/T)
    # Get limits for xi
    xi_min,xi_max = limits_xi(n_0,nus)
    xi_guess      = 0.5 * (xi_max + xi_min)
    result        = root_scalar(lambda xi: get_Qp_PT(n_0,nus,xi,P)-Kp_T,x0=xi_guess,bracket=(xi_min,xi_max),xtol=1E-8)
    xi_eq         = float(result.root)
    return xi_eq
# --------------------------------------------- #
def get_xieq_VT(V,T,n_0,nus,refdata):
    '''Extent of reaction in terms of V and T'''
    # Get Eq constant
    dGo_T = get_DGo(T,refdata)
    Kp_T  = np.exp(-dGo_T/R/T)
    # Get limits for xi
    xi_min,xi_max = limits_xi(n_0,nus)
    xi_guess      = 0.5 * (xi_max + xi_min)
    result        = root_scalar(lambda xi: get_Qp_VT(n_0,nus,xi,V,T)-Kp_T,x0=xi_guess,bracket=(xi_min,xi_max))
    xi_eq         = float(result.root)
    return xi_eq
# ============================================= #


# ============================================= #
# Functions for PART 2: STATIST.-THERMODYNAMICS #
# ============================================= #
def level_to_string(functional,basis):
    level      = rf"{functional.upper():s}_{basis.upper():s}"
    # replace * by _ast_ to avoid name problems
    level      = level.replace("**","_ast_ast_")
    level      = level.replace("*","_ast_")
    return level
# --------------------------------------------- #
def files_of_interest(molecule,functional="",basis="",DFTGRID=""):

    level      = level_to_string(functional,basis)
    # filenames
    try   : sgrid = rf".grid{DFTGRID:d}"
    except: sgrid = rf""
    xyz_guess  = rf"xyz_guess-{molecule:s}.xyz"
    xyz_opt    = rf"xyz_optim-{molecule:s}.{level:s}{sgrid:s}.xyz"
    output_opt = rf"pyscf_opt-{molecule:s}.{level:s}{sgrid:s}.out"
    output_frq = rf"pyscf_frq-{molecule:s}.{level:s}{sgrid:s}.out"
    # return filenames
    return xyz_guess,xyz_opt,output_opt,output_frq
# --------------------------------------------- #
def pubchem_cid(cid):
    '''PUBCHEM: get geom'''
    symbols,coords,smiles = None,None,None

    try   : compound = pcp.Compound.from_cid(int(cid))
    except: compound = None

    if compound is not None:
       print(rf"     - geometry retrieved!")
       print(rf"")
       print(rf"     - information:")
       print(rf"       * PubChem CID       : {compound.cid}")
       print(rf"       * molecular formula : {compound.molecular_formula:s}")
       # print(rf"       * charge            : {compound.charge:d}")
       print(rf"       * SMILES            : '{compound.smiles:s}'")
       print("")
       smiles  = compound.smiles
       # for diatomic molecules, 3D geoemtry is not stored in PubChem
       try   : geom = pcp.get_compounds(int(cid),"cid",record_type='3d')[0]
       except: geom = pcp.get_compounds(int(cid),"cid")[0]
       symbols = [atom.element           for atom in geom.atoms]
       # z-coordinate may be None (when record_type='3d' is not used)
       coords  = [tuple(0.0 if c is None else c for c in (atom.x, atom.y, atom.z)) for atom in geom.atoms]
    return symbols,coords,smiles
# --------------------------------------------- #
def rdkit_smiles2geom(smiles):
    '''RDKIT: get geom (from SMILES)'''

    symbols,coords = None,None
    try:
       m   = Chem.AddHs(Chem.MolFromSmiles(smiles))
       cid = AllChem.EmbedMolecule(m)
       if cid >= 0:
          symbols = [atom.GetSymbol() for atom in m.GetAtoms()]
          coords  = [list(m.GetConformer().GetAtomPosition(i)) for i in range(m.GetNumAtoms())]
          mformu  = {s:0 for s in symbols}
          for s in symbols: mformu[s] += 1
          mformu  = "".join([k+str(v) if v!=1 else k for k,v in mformu.items()])
          print(rf"     - geometry generated!")
          print(rf"")
          print(rf"     - information:")
          print(rf"       * molecular formula = {mformu:s}")
          print(rf"       * SMILES            = '{smiles:s}'")
          print("")
    except: pass

    return symbols,coords,smiles
# --------------------------------------------- #
def data_2_xyz(symbols,xcc,fname,smiles=""):
    '''data (from smiles) to .xyz file'''
    string  = rf"{len(symbols):.0f}"+"\n"
    string += rf"Cartesian coordinates for SMILES: {smiles:s}"+"\n"
    for idx,symbol in enumerate(symbols):
        xx,yy,zz = [coord for coord in xcc[idx]]
        string += rf"{symbol:2s}     {xx:9.5f}  {yy:9.5f}  {zz:9.5f}" + "\n"
    with open(fname,'w') as asdf: asdf.write(string)
    # download file button
    btn = download_file(fname)
    display(btn)
    print("")
# --------------------------------------------- #
def read_xyz(filename):
    with open(filename) as f: lines = f.readlines()
    nat = int(lines[0])
    symbols = []
    coords = []
    for line in lines[2:2+nat]:
        parts = line.split()
        if len(parts) < 4: continue
        sym, x, y, z = parts[:4]
        symbols.append(sym)
        coords.append([float(x), float(y), float(z)])
    return symbols, np.array(coords)
# --------------------------------------------- #
def show_indented(view, indent_px=40):
    display(HTML(f"""
        <div style="margin-left: {indent_px}px;">
            {view._make_html()}
        </div>
    """))
# --------------------------------------------- #
def create_visualization_xyz(xyz_file):
    # draw molecule
    view = py3Dmol.view(width=300, height=200)
    view.addModel(open(xyz_file, 'r').read(), 'xyz')
    view.setStyle({'stick': {'singleBonds': False}, 'sphere': {'scale': 0.3}})
    # Automatic labels: 0-based
    index_style = {"fontColor": "black","fontSize": 12,"showBackground": True ,"backgroundColor": "white","backgroundOpacity": 0.7}
    #elemn_style = {"fontColor": "black","fontSize": 12,"showBackground": False,"inFront": True,"alignment": "center","offset": {"x": 0, "y": 10, "z": 0}}
    view.addPropertyLabels("index",{},index_style)
    #view.addPropertyLabels("elem" ,{},elemn_style)
    #zoom
    view.zoomTo()
    view.zoom(2.5)
    return view
# --------------------------------------------- #
def pyscf_carryout_opt(molecule,unpaired,charge,functional,basis,DFTGRID,bsym=False):

    # Files of interest
    xyz_guess,xyz_opt,output_opt,output_frq = files_of_interest(molecule,functional,basis,DFTGRID)

    # Generate molecule from .xyz file
    with io.capture_output() as captured:
        mol = pyscf.gto.Mole()
        mol.atom         = xyz_guess
        mol.spin         = unpaired
        mol.charge       = charge
        mol.basis        = basis
        mol.output       = output_opt
        mol.verbose      = 4
        if bsym:
           mol.symmetry     = True
           mol.symmetrize   = True
           mol.symmetry_tol = 1e-2
        mol.build()

    # Define DFT method and mesh grid (from 0 to 9, default = 3)
    if unpaired == 0: mf = mol.RKS(xc=functional)
    else            : mf = mol.UKS(xc=functional)
    mf.grids.level = DFTGRID
    mf.max_cycle   = 200
    mf.conv_tol    = 1e-7

    # run SCF and avoid printing information
    print("     - single point calculation of guess geometry...")
    t1   = time.time()
    with io.capture_output() as captured: mf   = mf.run()
    Etot = mf.e_tot
    t2   = time.time()
    print(rf"       Etot(guess geometry) = {Etot:.5f} hartree")
    print(rf"       SCF took {t2-t1:.1f} seconds")

    # optimization [avoid printing information]
    t1   = time.time()
    print("     - geometry optimization...")
    conv_params = {}
    conv_params['convergence_energy'] = 5.0e-7  # Eh
    conv_params['convergence_gmax'  ] = 2.0e-4  # Eh/Bohr
    with io.capture_output() as captured:
       try:
          opt_geom = geometric_solver.optimize(mf, maxsteps=300, **conv_params)
       except:
          conv_params['convergence_energy'] = 1.0e-6  # Eh
          conv_params['convergence_gmax'  ] = 1.0e-4  # Eh/Bohr
          opt_geom = geometric_solver.optimize(mf, maxsteps=300, **conv_params)
    for line in opt_geom.tostring().split("\n"): print(7*" "+line)
    t2   = time.time()
    print(rf"       geometry opt. took {t2-t1:.1f} seconds")

    pyscf.gto.tofile(opt_geom,xyz_opt)
# --------------------------------------------- #
def pyscf_carryout_frq(molecule,unpaired,charge,functional,basis,DFTGRID,bsym=False):

    # Files of interest
    xyz_guess,xyz_opt,output_opt,output_frq = files_of_interest(molecule,functional,basis,DFTGRID)

    # Generate molecule from optimized geometry stored in xyz file
    with io.capture_output() as captured:
        mol = pyscf.gto.Mole()
        mol.atom         = xyz_opt
        mol.spin         = unpaired
        mol.charge       = charge
        mol.basis        = basis
        mol.output       = output_frq
        mol.verbose      = 6
        if bsym:
           mol.symmetry     = True
           mol.symmetrize   = True
           mol.symmetry_tol = 1e-2
        mol.build()

    # Define DFT method and mesh grid
    if unpaired == 0: mf = mol.RKS(xc=functional)
    else            : mf = mol.UKS(xc=functional)
    mf.grids.level = DFTGRID
    mf.max_cycle   = 200
    mf.conv_tol    = 1e-7

    print("     - Hessian calculation...")
    t1   = time.time()

    # run SCF and avoid printing information
    with io.capture_output() as captured: mf   = mf.run()
    Etot = mf.e_tot
    print(rf"       Etot(optim geometry) = {Etot:.5f} hartree")

    # Carry out Hessian calculation
    with io.capture_output() as captured: hessian = mf.Hessian().kernel()

    # add Hessian matrix to output_frq
    H4 = np.asarray(hessian)
    with open(output_frq, "a") as f:
        f.write("\n*** HESSIAN BY 3x3 ATOMIC BLOCKS (Hartree/Bohr^2) ***\n")
        for i in range(mol.natm):
            for j in range(mol.natm):
                f.write(f"\n# Block ({i+1},{j+1})  [atom {i+1} vs atom {j+1}]\n")
                np.savetxt(f, H4[i, j], fmt=" % .6e")

    t2   = time.time()
    print(rf"       Hessian calc. took {t2-t1:.1f} seconds")

    # return data
    return mol, mf, hessian
# --------------------------------------------- #
def pyscf_extract(mol, mf, hessian, unpaired):

    # translational info
    masses  = mol.atom_mass_list(isotope_avg=True)
    mass    = sum(masses)

    # vibrational info
    freqs       = thermo.harmonic_analysis(mf.mol, hessian)['freq_au']
    au2hz       = (1/2/np.pi)*(Eh/m_u/a_0**2)**0.5
    freqs_Hz    = [f*au2hz for f in freqs]
    wavenum_m   = [f/c_0   for f in freqs_Hz]

    info_thermo = thermo.thermo(mf, freqs) # by default, at 298.15 K and 101325 Pa
    ZPE         = info_thermo['ZPE'][0]

    # rotational info
    A,B,C = thermo.rotation_const(masses,mol.atom_coords(),unit='GHz')
    # if linear, two are equal, the other is infinity
    if   np.isinf(A) and abs(B-C) < ZERO3: A,B,C,linear = B,None,None,True
    elif np.isinf(B) and abs(A-C) < ZERO3: A,B,C,linear = C,None,None,True
    elif np.isinf(C) and abs(A-B) < ZERO3: A,B,C,linear = A,None,None,True
    else                                 :       linear =             False
    sigma = thermo.rotational_symmetry_number(mol)

    # electronic info
    E0 = info_thermo['E0' ][0]

    # collect data
    data = {}
    data["natoms"  ] = mol.natm
    data["mass"    ] = mass
    data["rotcons" ] = [i*1E9 if i is not None else i for i in [A,B,C]] # in Hz
    data["rotsigma"] = sigma
    data["islinear"] = linear
    data["freqs"   ] = wavenum_m # in 1/m
    data["ZPE"     ] = ZPE       # in hartree
    data["unpaired"] = unpaired
    data["E0"      ] = E0        # in hartree

    # return data
    return data
# --------------------------------------------- #
def optimize_and_freqs(molecule,unpaired,charge,functional,basis,DFTGRID,bsym=False):
    # files
    xyz_guess,xyz_opt,output_opt,output_frq = files_of_interest(molecule,functional,basis,DFTGRID)
    # check if calculation was already carried out
    args = (molecule,unpaired,charge,functional,basis,DFTGRID,bsym)
    # geometry optimization
    if not os.path.exists(xyz_opt): pyscf_carryout_opt(*args)
    # Hessian calculation
    mol, mf, hessian = pyscf_carryout_frq(*args)
    # Extract data
    return pyscf_extract(mol,mf,hessian,unpaired)
# --------------------------------------------- #
def pfn_translational(T,mass):

    beta    = 1/(k_B*T)

    # to SI
    mass_SI = mass * m_u

    # translational contribution
    V_per_molec = 1 / (beta * P_o)
    broglie_wvl = ((beta * h**2) / (2 * np.pi * mass_SI))**(0.5)
    q_translat  = V_per_molec / broglie_wvl**3

    # d(lnq_tr)/dbeta at constant volume
    dlnqdbeta_v = - 3/2 * (1/beta)

    # d(lnq_tr)/dbeta at constant pressure
    dlnqdbeta_p = - 5/2 * (1/beta)

    return q_translat, dlnqdbeta_v
# --------------------------------------------- #
def pfn_rotational(T,A,B,C,linear,sigma):

    beta = 1/(k_B*T)

    # rotational constants (Hz --> 1/m) and rot temperature
    A /= c_0; theta_A = (h*c_0/k_B) * A

    if linear:
      q_rotational = T / theta_A
      dlnqdbeta    = - (1/beta)
    else:
      B /= c_0; theta_B = (h*c_0/k_B) * B
      C /= c_0; theta_C = (h*c_0/k_B) * C
      q_rotational = (np.pi * T**3 / (theta_A * theta_B * theta_C))**(1/2)
      dlnqdbeta    = - 3/2 * (1/beta)

    q_rotational /= sigma

    return q_rotational, dlnqdbeta
# --------------------------------------------- #
def pfn_vibrational(T,freqs):

    beta = 1/(k_B*T)

    q_vibrational = 1.0
    dlnqdbeta     = 0.0
    for freq in freqs:
        nu    = freq * c_0 # in Hz
        theta = h*nu/k_B
        q_vibrational *= 1 / (1 - np.exp(-theta/T)) # from ZPE
        dlnqdbeta     += - h*nu / (np.exp(theta/T) - 1)

    return q_vibrational, dlnqdbeta
# --------------------------------------------- #
def pfn_electronic(T,unpaired):

    q_electronic = unpaired + 1
    dlnqdbeta    = 0.0

    return q_electronic, dlnqdbeta
# --------------------------------------------- #
def compute_thermodynamics(T,molecule,key,dftdata):

    # unpack data
    natoms   = dftdata[molecule][key]["natoms"  ]
    mass     = dftdata[molecule][key]["mass"    ]
    A,B,C    = dftdata[molecule][key]["rotcons" ]
    linear   = dftdata[molecule][key]["islinear"]
    sigma    = dftdata[molecule][key]["rotsigma"]
    freqs    = dftdata[molecule][key]["freqs"   ]
    unpaired = dftdata[molecule][key]["unpaired"]
    E0       = dftdata[molecule][key]["E0"      ]
    ZPE      = dftdata[molecule][key]["ZPE"     ]

    # Get partition functions
    q_translat   , dlnqtdbeta = pfn_translational(T,mass)
    if natoms == 1:
       q_rotational , dlnqrdbeta = 1.0, 0.0
       q_vibrational, dlnqvdbeta = 1.0, 0.0
    else:
       q_rotational , dlnqrdbeta = pfn_rotational(T,A,B,C,linear,sigma)
       q_vibrational, dlnqvdbeta = pfn_vibrational(T,freqs)
    q_electronic , dlnqedbeta = pfn_electronic(T,unpaired)
    q_tot         = q_translat * q_rotational * q_vibrational * q_electronic
    dlnqdbeta_tot = dlnqtdbeta + dlnqrdbeta + dlnqvdbeta + dlnqedbeta

    # Info line
    line  = rf" {q_translat:.4E} | {q_rotational:.4E} | {q_vibrational:.4E} | {q_electronic:.4E} | {q_tot:.4E} | {(E0+ZPE):12.5f}"

    # calculate U and S (in S.I.; per molecule)
    Eref = (E0+ZPE)*Eh
    U    = Eref - dlnqdbeta_tot
    S    = - dlnqdbeta_tot/T + k_B * np.log(q_tot * np.e)

    # calculate H and G (per molecule)
    H    = U + k_B * T
    G    = H - S   * T

    # return data in J and J/K
    return U, H, S, G, line
# --------------------------------------------- #
def sum_squared_errors(DCP,T,DGo_T,TREF,DHo_ref,DSo_ref):
    refdata   = (DHo_ref,DSo_ref,DCP*R,TREF)
    DGo_model = get_DGo(T,refdata)
    residuals = DGo_model - DGo_T
    return np.sum(residuals**2)
# --------------------------------------------- #
def calculate_rms(x1,x2):
    return (sum([(i-j)**2 for i,j in zip(x1,x2)])/len(x1))**0.5
# --------------------------------------------- #
def freq_to_nu_and_theta(freq_cm):
    nu    = 100*freq_cm*c_0
    theta = h*nu/k_B
    return nu, theta
# --------------------------------------------- #
def vib_contribution(freq_cm,T):
    nu,theta = freq_to_nu_and_theta(freq_cm)
    thetaT   = theta/T
    contri   = (thetaT * np.exp(-thetaT/2)/(1-np.exp(-thetaT))) **2
    return contri
# --------------------------------------------- #
def vib_contri_avera(freq_cm,T1,T2):
    nu,theta = freq_to_nu_and_theta(freq_cm)
    average  = 1/(np.exp(theta/T2)-1) - 1/(np.exp(theta/T1)-1)
    average  = average * theta/(T2-T1)
    return average
# ============================================= #


# ============================================= #
#      FUNCTIONS FOR PRINTING INFORMATION       #
# ============================================= #
def print_info_eq(magnitude,PVT_0,PVT_eq,molecules,xi_eq,n_eq,y_eq,p_eq,P_eq,Kp_v1,Kp_v2,Ky_v1,Ky_v2):

    P_0 ,V_0 ,T_0  = PVT_0
    P_eq,V_eq,T_eq = PVT_eq

    sP_0  = rf"{P_0/1E5:.2f}"
    sP_eq = rf"{P_eq/1E5:.2f}"
    Pformat = max(len(sP_0),len(sP_eq))

    sV_0  = rf"{V_0*1E3:.2f}"
    sV_eq = rf"{V_eq*1E3:.2f}"
    Vformat = max(len(sV_0),len(sV_eq))

    print("")
    print(fr"Initial  conditions: ({T_0:.2f}K,{sP_0:{Pformat:d}s}bar,{sV_0:{Vformat:d}s}L)")
    print("")
    print(fr"Equilib. conditions: ({T_eq:.2f}K,{sP_eq:{Pformat:d}s}bar,{sV_eq:{Vformat:d}s}L)")
    print("")
    print(fr"   ==> equilibrium found at xi_eq = {xi_eq:.4f} mol")
    print("")

    print(fr"   Number of moles & molar fraction at equilibrium (from xi_eq):")
    for j,molecule in enumerate(molecules):
        n_eq_j = n_eq[j]
        y_eq_j = y_eq[j]
        print(fr"   * n_i = {n_eq_j:6.4f} mol, y_i = {y_eq_j:7.4f} (i = {molecule:s})")
    print("")

    print(fr"   Partial pressures at equilibrium:")
    for j,molecule in enumerate(molecules):
        p_eq_j = p_eq[j]/1E5
        print(fr"   * p_i = {p_eq_j:7.4f} bar (i = {molecule:s})")
    print("")


    if 9999 > Kp_v1 > 0.100: sformat1 = ".3f"
    else                   : sformat1 = ".2E"
    if 9999 > Ky_v1 > 0.100: sformat2 = ".3f"
    else                   : sformat2 = ".2E"

    reldiff = 100*abs(Kp_v2-Kp_v1)/Kp_v1
    if reldiff < 5.0:
        print(fr"   Value of Kp:")
        print(fr"   * from Delta_r{{G}}^o --> Kp = {Kp_v2:{sformat1}} [*]")
        print(fr"   * from p_i values   --> Kp = {Kp_v1:{sformat1}}")
        print("")
        print(fr"   Value of Ky:")
        if Ky_v2 is not None:
          print(fr"   * from Delta_r{{G}}^* --> Ky = {Ky_v2:{sformat2}} [*]")
        print(fr"   * from y_i values   --> Ky = {Ky_v1:{sformat2}}")
        print("")
        print(fr"   The previous values for the equilibrium constants may vary slightly")
        print(fr"   due to numerical errors. Trust the values with the [*].")
    else:
        print(fr"   Value of Kp:")
        print(fr"   * from Delta_r{{G}}^o --> Kp = {Kp_v2:{sformat1}}")
        print("")
        if Ky_v2 is not None:
          print(fr"   Value of Ky:")
          print(fr"   * from Delta_r{{G}}^* --> Ky = {Ky_v2:{sformat2}}")
    print("")
# --------------------------------------------- #
def print_sym_nums(MOLECULES,LEVELS,DFTDATA):
    nn   = max([len(molecule) for molecule in MOLECULES])
    smol = "Molecule"
    nn   = max(len(smol),nn)
    while len(smol) < nn: smol = " "+smol+" " 
    if    len(smol) > nn: smol = smol[:-1]

    SLEVELS = [rf"{functional}/{basis}" for functional,basis in LEVELS]
    mm      = max([len(slevel) for slevel in SLEVELS])
    SLEVELS = [("%%%is"%mm)%slevel for slevel in SLEVELS]

    line = rf" {smol} | "+ " | ".join([slevel for slevel in SLEVELS])
    divi = "-"*len(line)
    print(line)
    print(divi)
    for molecule in MOLECULES:
        if molecule not in DFTDATA: continue
        line = rf" {molecule:{nn}s} "
        for functional,basis in LEVELS:
            key = (functional,basis)
            if key in DFTDATA[molecule]:
               sigma  = DFTDATA[molecule][key]["rotsigma"]
               ssigma = str(sigma)
               while len(ssigma) < mm: ssigma = " "+ssigma+" "
               if    len(ssigma) > mm: ssigma = ssigma[:-1]
               line += rf"| {ssigma} "
            else: line += "|  "+mm*" "
        print(line)
    print(divi)
    print("")
# --------------------------------------------- #
def geometric_info_xyz(xyz_file,geominfo):
    info = ""
    symbols, coords = read_xyz(xyz_file)
    for ii in geominfo:
        if len(ii) == 2:
          at1,at2 = ii
          distance = np.linalg.norm(coords[at1] - coords[at2])
          sbond    = rf"{symbols[at1]:s}{at1:d}-{symbols[at2]:s}{at2:d}"
          info    += rf"       * d({sbond:8s}) = {distance:.3f} Å" + "\n"
        elif len(ii) == 3:
          at1,at2,at3 = ii
          v1 = coords[at1] - coords[at2]
          v2 = coords[at3] - coords[at2]
          cosang = np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2))
          cosang = np.clip(cosang, -1.0, 1.0)
          angle  = np.degrees(np.arccos(cosang))
          sangle = rf"{symbols[at1]:s}{at1:d}-{symbols[at2]:s}{at2:d}-{symbols[at3]:s}{at3:d}"
          info  += rf"       * ∠({sangle:8s}) = {angle:.1f}°" + "\n"
    return info
# --------------------------------------------- #
def pyscf_printdata(data):
    '''print information of interest after Hessian calc'''
    # unpack data
    mass     = data["mass"    ]
    A,B,C    = data["rotcons" ]
    linear   = data["islinear"]
    sigma    = data["rotsigma"]
    freqs    = data["freqs"   ]
    ZPE      = data["ZPE"     ]
    unpaired = data["unpaired"]
    E0       = data["E0"      ]

    # print fata
    INFO  = rf"     - total mass (amu)              : {mass:.2f}" + "\n"
    if linear: INFO += rf"     - rotational constant  (GHz)    : {(A*1E-9):.2f}" + "\n"
    else     : INFO += rf"     - rotational constants (GHz)    : " + "  ".join(["%.2f"%(ii*1E-9) for ii in (A,B,C)]) + "\n"
    INFO += rf"     - rotational symmetry number    : {sigma:d}" + "\n"
    INFO += rf"     - vibrational wavenumbers (1/cm):"+"\n"
    # frequencies are, actually, wavenumbers in (1/m)
    for i in range(0,len(freqs),5):
        INFO += "         "+"  ".join(["%8.2f"%(ii/100) for ii in freqs[i:i+5]])+"\n"
    INFO += rf"     - zero point energy (hartree)   : {ZPE:.5f}"+"\n"
    print(INFO)
# ============================================= #


# ============================================= #
# ----          PLOTTING FUNCTIONS         ---- #
# ============================================= #
def plot_DGo_T(T,TREF,refdata):
    # Calculation of DG at Tref
    DGo_ref = get_DGo(TREF,refdata)
    Keq_ref = np.exp(-DGo_ref/R/TREF)

    # Calculation of DG at other temperatures
    DGo_T = get_DGo(T,refdata)
    Keq_T = np.exp(-DGo_T/R/T)

    # Plot results
    plt.rcParams['text.usetex'] = True

    # Create figure with two panels side by side
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    # LEFT PANEL #
    labeldot = r"$\Delta_{r}{G}^\circ(T_{\rm ref}) = %.2f \;\; \mathrm{kJ/mol}$"%(DGo_ref/1000)
    ax[0].plot(T   ,DGo_T  /1000,'k-',zorder=1)
    ax[0].plot(TREF,DGo_ref/1000,'ro',zorder=3,label=labeldot)
    ax[0].tick_params(axis='both', labelsize=FONTSIZE[4])
    ax[0].set_xlabel(r'$T \;\; (\mathrm{K})$'                        ,fontsize=FONTSIZE[5])
    ax[0].set_ylabel(r'$\Delta_{r} G^{\circ} \;\; (\mathrm{kJ/mol})$',fontsize=FONTSIZE[5])
    ax[0].legend(loc="best",fontsize=FONTSIZE[3])

    # RIGHT PANEL #
    if 1E-2 < Keq_ref < 1000: sKeq_ref = "%.3f"%Keq_ref
    else                    : sKeq_ref = "%.3e"%Keq_ref
    labeldot = r"$K_{p}^\circ(T_{\rm ref}) = %s$"%sKeq_ref
    ax[1].plot(T   ,Keq_T  ,'k-',zorder=1)
    ax[1].plot(TREF,Keq_ref,'ro',zorder=3,label=labeldot)
    ax[1].tick_params(axis='both', labelsize=FONTSIZE[4])
    ax[1].set_xlabel(r'$T \;\; (\mathrm{K})$'     ,fontsize=FONTSIZE[5])
    ax[1].set_ylabel(r'$K_{p}^\circ(T_{\rm ref})$',fontsize=FONTSIZE[5])
    ax[1].legend(loc="best",fontsize=FONTSIZE[3])

    # ---- download button ----
    fname = 'plot_DGo_T.svg'
    fig   = plt.gcf()
    fig.savefig(fname, bbox_inches='tight')
    btn   = download_file(fname)

    # ---- show plot and close ----
    display(btn)
    plt.show()
    plt.close()

    return DGo_T
# --------------------------------------------- #
def plot_DGo_T_statmech(T,DGo_T,DGo_model,DCP,TREF,DGo_ref,key):

    plt.rcParams['text.usetex'] = True

    if DCP is not None: DCP = DCP/R
    # Plot results
    plt.plot(T,[ii/1000 for ii in DGo_T]    ,'k-',zorder=1,label=r"Stat-Mech")

    if DGo_model is not None:
       plt.plot(T,[ii/1000 for ii in DGo_model],'r--',zorder=1,label=r"Class-Therm with $\Delta_{r}C_{p}^\circ= %.2f\cdot R$"%DCP)

    if TREF is not None:
       plt.plot(TREF,DGo_ref/1000,'ro')

    # Format plot
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(-15,25)

    plt.xlabel(r'$T \;\; (\mathrm{K})$'                        ,fontsize=16)
    plt.ylabel(r'$\Delta_{r} G^{\circ} \;\; (\mathrm{kJ/mol})$',fontsize=16)

    plt.legend(loc="best",fontsize=12)

    if DGo_model is not None:
       level = level_to_string(key[0],key[1])
       fname = rf"plot_DGo_T_stat_{level:s}_{DCP:5.2f}.svg"
       # ---- download button ----
       fig = plt.gcf()
       fig.savefig(fname, bbox_inches='tight')
       btn = download_file(fname)
       display(btn)

    # ---- show plot and close ----
    plt.show()
    plt.close()
# --------------------------------------------- #
def plot_gibbshelmholtz(T,DGo_T,refdata):
    # --- DH^o at each temperature ---
    DHo_T = get_DHo(T,refdata)
    yy    = -DHo_T / (R * T**2)

    # Plot results
    fig, ax = plt.subplots()
    ax.plot(T,yy,'k-',zorder=1,label=r'$(a) \;\; f=-\Delta_r H^{\circ}/T^2$ \quad\quad [analytic]')

    # Format plot
    ax.tick_params(axis='x', labelsize=FONTSIZE[4])
    ax.tick_params(axis='y', labelsize=FONTSIZE[4])

    ax.set_xlabel(r'$T \;\; (\mathrm{K})$'       , fontsize=FONTSIZE[5])
    ax.set_ylabel(r'$f/R \;\; (\mathrm{K}^{-1})$', fontsize=FONTSIZE[5])

    # --- Calculate numerical derivative ---
    numslope = np.gradient(DGo_T / (R*T) , T)

    # --- Add data to plot (removing first and last points, due to numerical error) ---
    ax.plot(T[1:-1],numslope[1:-1],'xr',zorder=2, markeredgewidth=2.5,label=r'$(b)  \;\; f = d/dT(\Delta_r G^{\circ}/T)$ \;\, [numeric]')

    ax.legend(loc="best",fontsize=FONTSIZE[3])

    # ---- download button ----
    fname = "plot_gibbshelmholtz.svg"
    fig   = plt.gcf()
    fig.savefig(fname, bbox_inches='tight')
    btn   = download_file(fname)

    # Show plot, button and close
    display(btn)
    display(fig)
    plt.close()
# --------------------------------------------- #
def plot_DG_PT(T,P, fixed_args):

    molecules,nus,n_0,xis,refdata = fixed_args

    V_0 = n_0.sum()*R*T/P

    # ---- Calculate DG* ----
    DGo   = get_DGo(T,refdata)
    Gast  = get_Gast(T,P,nus,refdata)

    # ---- Terms in DGtot ----
    DGlin  = Gast * xis
    DDGmix = get_DDGmix(xis,T,P,n_0,nus)
    DGtot  = DGlin + DDGmix

    # ---- minimum of DGtot ---
    xi_eq  = get_xieq_PT(P,T,n_0,nus,refdata)
    minDG  = get_G_PT(xi_eq,P,T,n_0,nus,refdata)
    minDG  = minDG/(R*T)

    # ---- minimum of DGmix ---
    min2xi  = xis[np.argmin(DDGmix)]
    min2DG  = np.min(DDGmix)/(R*T)

    # ---- Number of moles and pressure at equilibrium ----
    n_eq = n_0 + xi_eq * nus
    y_eq = n_eq / n_eq.sum()
    p_eq = P * y_eq
    V_eq = (n_eq.sum())*R*T/P

    # ---- Calculation of Kp and Ky ----
    # (a) using molar fractions and partial pressures (from xi_eq)
    Kp_v1, Ky_v1 = 1.0, 1.0
    for p_eq_j,nu_j in zip(p_eq,nus): Kp_v1 *= (p_eq_j/P_o)**nu_j
    for y_eq_j,nu_j in zip(y_eq,nus): Ky_v1 *= y_eq_j**nu_j
    # (b) using deltaG values
    Ky_v2 = np.exp(-Gast/(R*T))
    Kp_v2 = np.exp(-DGo /(R*T))
    # (c) directly from xi_eq (specific case A --> 2B)
    Ky_v3 = 4*xi_eq * xi_eq / (1-xi_eq*xi_eq)

    # --- Plot data ---
    fig = plt.figure(figsize=(5.50,3.67))
    plt.plot(xis, DGlin/(R*T), color='b',ls="--",label=r"$f = \xi \cdot \Delta_{r} G^{\ast}$")
    plt.plot(xis,DDGmix/(R*T), color='r',ls=":" ,label=r"$f = \Delta_{\rm mix}G(\xi) - \Delta_{\rm mix}G(0)$")
    plt.plot(xis, DGtot/(R*T), color='k',ls="-" ,label=r"$f = G(\xi) - G(0)$")
    plt.plot( xi_eq ,minDG ,'ko' )
    plt.plot( min2xi,min2DG,'rx' )

    # some formatting
    plt.xticks(fontsize=FONTSIZE[1])
    plt.yticks(fontsize=FONTSIZE[1])
    plt.xlabel(r'$\xi \;\; \mathrm{[mol]}$'           , fontsize=FONTSIZE[2])
    plt.ylabel(r'$f \; / \; (RT) \;\; \mathrm{[mol]}$', fontsize=FONTSIZE[2])
    plt.title(fr'T={T:.0f} K, p={P/1E5:.2f} bar')

    # secondary grid in x-axis
    ax = plt.gca()
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.1))
    ax.grid(True, which='minor', axis='x', alpha=0.15)
    ax.tick_params(axis='x', which='minor',bottom=False,top=False,length=0)
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.legend(loc="best", fontsize=FONTSIZE[0])

    # --- update global variable: last_fig ---
    global last_fig
    last_fig = plt.gcf()

    # --- Show and close figure ---
    plt.show()
    plt.close()

    # ---- Print info ----
    PVT_0  = (P,V_0 ,T)
    PVT_eq = (P,V_eq,T)
    # print(fr"Value for Delta_r{{G}}^o({T:.2f} K) = {DGo/1000:.2f} kJ/mol")
    # print(fr"Value for Delta_r{{G}}^*({T:.2f} K) = {Gast/1000:.2f} kJ/mol")
    print_info_eq("G",PVT_0,PVT_eq,molecules,xi_eq,n_eq,y_eq,p_eq,P,Kp_v1,Kp_v2,Ky_v1,Ky_v2)
# --------------------------------------------- #
def plot_DA_VT(T,V, fixed_args):

    molecules,nus,n_0,xis,refdata = fixed_args

    # ---- Terms in DAtot ----
    DAtot, term_lineal, term_nolin, P_xi = get_A_VT(xis,V,T,n_0,nus,refdata)
    # ---- minimum of DAtot ---
    idx_eq = np.argmin(DAtot)
    xi_eq  = xis[idx_eq]
    P_eq   = P_xi[idx_eq]
    minDA  = np.min(DAtot)/(R*T)

    # ---- minimum of DAtot ---
    xi_eq  = get_xieq_VT(V,T,n_0,nus,refdata)
    minDA  = get_A_VT(xi_eq,V,T,n_0,nus,refdata)[0]/(R*T)

    # ---- minimum of non-lineal term ---
    min2xi  = xis[np.argmin(term_nolin)]
    min2yy  = np.min(term_nolin)/(R*T)

    # ---- Number of moles and pressure at equilibrium ----
    n_eq = n_0 + xi_eq * nus
    y_eq = n_eq / n_eq.sum()
    p_eq = P_eq * y_eq

    # ---- Calculation of Kp and Ky ----
    # (a) using molar fractions and partial pressures (from xi_eq)
    Kp_v1, Ky_v1 = 1.0, 1.0
    for p_eq_j,nu_j in zip(p_eq,nus): Kp_v1 *= (p_eq_j/P_o)**nu_j
    for y_eq_j,nu_j in zip(y_eq,nus): Ky_v1 *= y_eq_j**nu_j
    # (b) using deltaG values
    DGo          = get_DGo(T,refdata)
    Kp_v2, Ky_v2 = np.exp(-DGo /(R*T)), None

    # --- Plot data ---
    plt.figure(figsize=(5.50,3.67))

    plt.plot(xis,term_lineal/(R*T), color='b',ls="--",label=r"$f = \xi \cdot \Delta_{r} G^{o}$")
    plt.plot(xis,term_nolin /(R*T), color='r',ls=":",label=r"$f = V \sum_i \left[ p_i \ln\frac{p_i}{e p^\circ} - p_i(0) \ln \frac{p_i(0)}{e p^\circ} \right]$")
    plt.plot(xis,      DAtot/(R*T), color='k',ls="-" ,label=r"$f = A(\xi) - A(0)$")

    plt.plot( xi_eq ,minDA ,'ko' )
    plt.plot(min2xi ,min2yy ,'rx' )

    plt.xticks(fontsize=FONTSIZE[1])
    plt.yticks(fontsize=FONTSIZE[1])
    plt.xlabel(r'$\xi \;\; \mathrm{[mol]}$'           , fontsize=FONTSIZE[2])
    plt.ylabel(r'$f \; / \; (RT) \;\; \mathrm{[mol]}$', fontsize=FONTSIZE[2])
    plt.title(fr'T={T:.0f} K,  V={V*1000:.2f} L')

    # secondary grid in x-axis
    ax = plt.gca()
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.1))
    ax.grid(True, which='minor', axis='x', alpha=0.15)
    ax.tick_params(axis='x', which='minor',bottom=False,top=False,length=0)
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.legend(loc="best",fontsize=FONTSIZE[0])

    # --- update global variable: last_fig ---
    global last_fig
    last_fig = plt.gcf()

    # --- Show and close figure ---
    plt.show()
    plt.close()

    # ---- Print info ----
    P_0 = n_0.sum() * R*T/V
    PVT_0  = (P_0 ,V,T)
    PVT_eq = (P_eq,V,T)
    print_info_eq("A",PVT_0,PVT_eq,molecules,xi_eq,n_eq,y_eq,p_eq,P_eq,Kp_v1,Kp_v2,Ky_v1,Ky_v2)
# --------------------------------------------- #
def plot_vib_average(Tmin,Tmax,freqmin,freqmax,freq=None):
    freqs    = np.linspace(freqmin,freqmax,51)
    averages = vib_contri_avera(freqs,Tmin,Tmax)
    yy_inf   = vib_contribution(freqs,Tmin)
    yy_sup   = vib_contribution(freqs,Tmax)
    plt.plot(freqs,yy_inf  ,'r--',label=rf"$T$={Tmin:.2f} K",zorder=2)
    plt.plot(freqs,yy_sup  ,'b--',label=rf"$T$={Tmax:.2f} K")
    plt.plot(freqs,averages,'k-' ,label=rf"average")
    plt.xlabel(r"Vib. frequency (cm$^{-1}$)"     ,fontsize=14)
    plt.ylabel(r"Vib. contribution, $n^{V^\ast}$",fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # --- update global variable: last_fig ---
    global last_fig
    last_fig = plt.gcf()

    if freq is not None:
       # limits so far
       xlim = plt.gca().get_xlim()
       ylim = plt.gca().get_ylim()
       # get average contribution for selected frequency
       aver =  vib_contri_avera(freq,Tmin,Tmax)
       # plot data for selected frequency
       plt.plot([freq,freq],[ylim[0],aver],'--',color="grey",zorder=1)
       plt.plot([xlim[0],freq],[aver,aver],'--',color="grey",zorder=1)
       plt.plot(freq,aver,'o',color="grey",label=rf"$n^{{V^\ast}} = {aver:.3f}$")
       # keep origonal limits
       plt.xlim(xlim)
       plt.ylim(ylim)
    plt.legend(loc="best",fontsize=11)

    # --- Show and close figure ---
    plt.show()
    plt.close()
# ============================================= #
