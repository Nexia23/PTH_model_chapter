import seaborn as sns
import numpy as np
import pandas as pd
import tellurium as te 
import re
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc
from matplotlib.cm import ScalarMappable
from colorsys import rgb_to_hsv,hsv_to_rgb
from collections import OrderedDict
import warnings
from functools import reduce
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')
import sys
sys.path.append('Estimation/')
from estimation import set_model_to_ss

def simulate_model(model_path: str = 'LCT_model/general_LCT_OIE.ant', 
                   infection_pars: dict = {}, name: str='general',
                   ss_duration:float=4.,pre_t:float=10.0, simulation_end:float=100.0, 
                   bool_med= True, bool_set_pars=False,
                   selections:list=['time', '[R]', '[iE]', '[E]','[M]','[P]', 'oiE', 'Hkt',
                                     'Hb', 'LDH', 'parasitemia', 'RPI', 'J_oiE_death',
                                     'J_P_birth', 'J_P_death',
                                     'k_E_death','rpi_step_func'] + [f'[oiE_{i}]' for i in range(1, 7)]):
    """Simulates model for given BH_max value and returns dataframe with results. 
    Includes presimulation for specific time (infection_after)."""
    model = te.loada(model_path)
   
    if bool_set_pars:
        for p in infection_pars:
            model.setValue(p, infection_pars[p])
    # ss simulation
    model = set_model_to_ss(model, infection_pars, model_name=name)
    
    model.M = 0
    model.events_medication_on = False 
    ss_res = model.simulate(-ss_duration-pre_t, -pre_t, 10, selections=selections)
    ss_res_df = pd.DataFrame(ss_res, columns=ss_res.colnames)    

    # infection simulation
    if 'M' in infection_pars.keys():
        model.M = infection_pars['M']
    else: 
        model.M = 1e3

    for p in infection_pars:
        model.setValue(p, infection_pars[p])
    
    # ACT addition
    model.events_medication_on = bool_med
    act_res = model.simulate(-pre_t, simulation_end, (pre_t+simulation_end+1), selections=selections) #int(simulation_end+pre_t)*10
    act_res_df = pd.DataFrame(act_res, columns=act_res.colnames) 

    # combine results
    res_df = pd.concat([ss_res_df, act_res_df], axis=0)    
    return res_df 

def t_relax(t,x):

    eq = get_steady_state(t,x)
    pos,h = find_peaks(x,eq)
    h = h["peak_heights"]
    A = np.abs((h - eq))
    A_i = np.argwhere(A < 0.01)
    
    if np.isfinite(eq) and A_i.shape[0] > 0:
        return t[pos[A_i[0,0]]]
    else:
        return np.inf
    

def get_steady_state(t,x):

    x_e = x[-int(0.1 * len(t)):]
    pos, h = find_peaks(x_e,np.mean(x_e) * 1.01)
    if len(pos) > 0:
        return np.inf
    else:
        return np.mean(x_e)

def max_rel_amplitude(t,x):

    eq = get_steady_state(t,x)
    if eq<=1e-6: eq=np.mean(x)
    pos,h = find_peaks(x,eq)
    if len(pos) == 0:
        return None,0
    i = np.argmax(h["peak_heights"])
    return pos[i], (h["peak_heights"][i] - eq)/eq


def max_abs_amplitude(t,x):
    n = 20 #rolling period
    local_max_vals = x.loc[x == x.rolling(n, center=True).max()]
    #df.loc[df[species] == df[species].rolling(n, center=True).min()]
    #print(local_max_vals)
    eq = get_steady_state(t,x)
    pos,h = find_peaks(x,eq)
    #print(pos, h)
    if len(pos) == 0:
        return None,0
    i = np.argmax(h["peak_heights"])
    return local_max_vals.max()#pos[i], h["peak_heights"][i]


def max_peak_width(t,x):
    #df.loc[df[species] == df[species].rolling(n, center=True).min()]
    #print(local_max_vals)
    eq = get_steady_state(t,x)
    pos,h = find_peaks(x,eq)

    #print(pos, h)
    if len(pos) == 0:
        return np.nan
    i = np.argmax(h["peak_heights"])
    local_max_vals = x.loc[0.95*h["peak_heights"][i]<x]
    return len(local_max_vals)#pos[i], h["peak_heights"][i]


def min_peak_width(t,x):
    n = 20 #rolling period
    x = x.reset_index(drop=True)
    local_min_vals = x.loc[x == x.rolling(n, center=True).min()]
    #df.loc[df[species] == df[species].rolling(n, center=True).min()]
    #print(local_min_vals)
    if len(local_min_vals)==0:
        return np.nan
    if local_min_vals.max()<=1e-6: return np.nan
    i = local_min_vals.idxmax()

    threshold = 1.05*local_min_vals.max()
    filtered_values =[]
    # Traverse left from the minimum value index
    left_index = i
    #print(x.iloc[left_index])
    while left_index >= 0 and x.iloc[left_index] <= threshold:
        filtered_values.append(left_index)
        left_index -= 1

    # Traverse right from the minimum value index
    right_index = i + 1
    while x.iloc[right_index] <= threshold:
        #print(x[right_index], threshold)
        filtered_values.append(right_index)
        right_index += 1
        if right_index == 510:break

    return len(filtered_values)#pos[i], h["peak_heights"][i]


def min_abs_amplitude(t,x):
    n = 20 #rolling period
    local_min_vals = x.loc[x == x.rolling(n, center=True).min()]
    if len(local_min_vals)==0:
        return np.nan
    if local_min_vals.max()<=1e-6: return [],None
    i = local_min_vals.idxmin()
    if len(local_min_vals) == 0:
        return None,0
    
    return i, local_min_vals.min()


def min_rel_amplitude(t,x):
    n = 20 #rolling period
    local_min_vals = x.loc[x == x.rolling(n, center=True).min()]
    if len(local_min_vals)==0:
        return np.nan
    eq = get_steady_state(t,x)
    if local_min_vals.max()<=1e-6: return np.nan
    i = local_min_vals.idxmin()
    if len(local_min_vals) == 0:
        return None,0
    print(local_min_vals)
    return i, (local_min_vals.min() -eq)/eq

    
def freq(t,x):

    eq = np.mean(x)
    pos,h = find_peaks(x,eq * 1.1)
    return 1/np.mean(np.diff(t[pos]))

def t_relax_smooth(t,x):
    eq = get_steady_state(t,x)
    def g(x,a,tau):
        return a * np.exp(-1/tau*x) + eq

    pos,h = find_peaks(x,eq)
    if len(pos) <= 1:
        return np.inf
    popt,pcov = curve_fit(g,x[pos],h["peak_heights"])
    if pcov[1,1]**2 > 1:
        return np.inf
        
    return popt[1]


def my_find_peaks(X,n=100,eps = 0.01):

    max_peaks = []
    min_peaks = []
    
    def p(x):
        x = (x - np.min(x))
        pos,h = find_peaks(x,eps)
        if len(pos)>0:
            max_peaks.append((x.index[pos].astype(int).tolist()))
            
        x = -x
        x = (x - np.min(x))
        pos,h = find_peaks(x, eps,)
        if len(pos)>0:
            min_peaks.append((x.index[pos].astype(int).tolist()))
        return 1
    
    X.rolling(n, center=True).apply(p)
    try:
        max_peaks = reduce(lambda a,b: a+b, max_peaks)
        max_peaks = np.unique(max_peaks)
    except TypeError:pass

    try:
        min_peaks = reduce(lambda a,b: a+b, min_peaks)
        min_peaks = np.unique(min_peaks)    
    except TypeError:pass

 
    return min_peaks,max_peaks


def peak_width(x, loc:list=[], type:str='max'):
 
    if len(loc)==0:
        return np.nan
    loc = int(loc[0])

    if x[loc]<=1e-6: return np.nan
    
    threshold = 1.01*x[loc]
    if type=='max':
        x =-x
        threshold = .95*x[loc] 
    
    filtered_values =[]
    # Traverse left from the minimum value index
    left_index = loc
    while left_index >= 0 and x[left_index] <= threshold:
        filtered_values.append(left_index)
        left_index -= 1

    # Traverse right from the minimum value index
    right_index = loc + 1
    while right_index< len(x) and x[right_index] <= threshold:
        #print(x[right_index], threshold)
        filtered_values.append(right_index)
        right_index += 1

    return len(filtered_values)

def new_max_peaks(t,x):
    
    x.index = pd.RangeIndex(len(x))    
    t.index = pd.RangeIndex(len(x))    

    min_peaks,max_peaks = my_find_peaks(x)

    result = {
        "maximum_magnitude":{f"peak_{i}":float(x.loc[loc]) for i,loc in enumerate(max_peaks)},
        "maximum_timepoint":{f"peak_{i}":float(t.loc[loc]) for i,loc in enumerate(max_peaks)},
        "maximum_width":{f"peak_{i}":float(peak_width(x,[loc],type='max')) for i,loc in enumerate(max_peaks)},
        "maximum_n_peaks":len(max_peaks),

        "minimum_magnitude":{f"peak_{i}":float(x.loc[loc]) for i,loc in enumerate(min_peaks)},
        "minimum_timepoint":{f"peak_{i}":float(t.loc[loc]) for i,loc in enumerate(min_peaks)},
        "minimum_width":{f"peak_{i}":float(peak_width(x,[loc],type='min')) for i,loc in enumerate(min_peaks)},
        "minimum_n_peaks":len(min_peaks),
    }
    return result


readouts = {
    #"frequency":freq,# mean frequency of oscillations that exceed the mean by 10%
    #"mean":lambda t,x:np.mean(x), # mean value
    'maximum_height': new_max_peaks,
    #'max_peak_width': lambda t,x: max_peak_width(t,x),
    #'min_peak_width': lambda t,x: min_peak_width(t,x),
    #"steady_state":lambda t,x: get_steady_state(t,x), # steady state value if solutin is non oscillatory
    #"max_rel_amplitude":lambda t,x: max_rel_amplitude(t,x)[1], # maximum peak amplitude as a multiple of the steady state
    #"max_abs_amplitude":lambda t,x: max_abs_amplitude(t,x), # max peak amplitude from 0

    #"min_rel_amplitude":lambda t,x: min_rel_amplitude(t,x), # maximum peak amplitude as a multiple of the steady state
    #"min_abs_amplitude":lambda t,x: min_abs_amplitude(t,x), # maximum peak amplitude as a multiple of the steady state

    #"t_relax":t_relax, # time until osciallations die out
    #"tau":lambda x,t: t_relax_smooth(x,t)# decay constant of the enveloping exponential function
}
#2024_02_19_14_16_1
update_pth= {
        "k_E_infect": 1.0351927708026337e-06,
        "tropism": 2.4824267001740767,
        "M": 100.0000740303447,
        "a_P_d": 9937191.134916836,
        "k_P_d": 0.009451056983160715,
        "r_P_d": 7.4265726002365495,
        "fac_R_d": 7.068386037799187e-10,
        "k_P_art_max": 2.709891804855283,
        "t_mat_P": 6.138028732895664,
        "k_iE_pit_frac": 0.6625125905514542,
        "t_E_death": 80.00000022727542,
        "s_BH": 7.105938524287848e-07, #change
        "LDH": 226.23256587878646,
        "k_M_death": 30.010891278063838,
        
        "Hkt_init":0.45
}
params_bounds = OrderedDict({
    's_BH': (5e-7, 5e-6),
    # 'k_E_infect': (1e-7, 2.5e-6),
    "Hkt_init":(0.35,0.55),
    #'slope_rpi':(1,200),      # own idea
    # 'a_P_d':(1e4,6e8),
    # 'k_M_death':(20,96),      # paper value 48 /day
    # 'k_iE_pit_frac': (0.1,1.1),
    # 'tropism': (0,100),
    # 't_mat_P':(2,14),        # medicine fact
    # 't_E_death':(70,130),

    #'r_P_d':(3,9.8),
    #'k_P_birth': (5e1, 1e3),
    
})
species_to_analyze = OrderedDict({
    'parasitemia' : ['parasitemia (%)', 1],
    #'[M]' : ['M (1e4/µl)', 1e4],
    '[iE]': ['iE (1e4/µl)', 1e4],
    'oiE': ['oiE (1e4/µl)', 1e4],
    '[E]': ['E (1e4/µl)', 1e4],
    '[R]': ['R (1e4/µl)', 1e4],
    '[P]': ['P (1e0/µl)', 1e0],
    #'J_P_birth':['J_P_birth',1e0],
    'Hb': ['Hb (g/dl)', 1],
    'Hkt': ['Hct', 1],
    'RPI': ['RPI', 1],
    'LDH': ['LDH (U/l)',1],
})
pre_t = 10
df = []

for p_name in params_bounds.keys(): 
    infection_dict = update_pth.copy()
    
    for m in np.logspace(-.1,.1,50):
        result = simulate_model(infection_pars= infection_dict | {p_name:update_pth[p_name]*m}, bool_med=True, bool_set_pars=True, 
                                pre_t=pre_t, simulation_end=80,selections=['time']+list(species_to_analyze.keys()))

        t = result['time']
        tester = result.div(result.iloc[0])
        tester = tester.drop(['time'],axis=1)
        tester[['parasitemia', '[iE]', 'oiE']] = result[['parasitemia', '[iE]', 'oiE']]      

        for species in result.columns[1:]:
            y = tester[species]
            for k,v in readouts.items():
                r = v(t,y)
                if isinstance(r,dict):
                    for sub_k,sub_r in r.items():
                        if isinstance(sub_r,dict):
                            for sub_sub_k,sub_sub_r in sub_r.items():
                                df.append({"readout_name":sub_k,"value":sub_sub_r,"sub_readout_name":sub_sub_k} | {"parameter":p_name, "parameter_value":m,"species":species})    
                        else:
                            df.append({"readout_name":k,"value":sub_r,"sub_readout_name":sub_k} | {"parameter":p_name, "parameter_value":m,"species":species})
                else:
                    try:
                        r = float(r)
                    except:continue
                    df.append({"readout_name":k,"value":r} | {"parameter":p_name, "parameter_value":m,"species":species})
df = pd.DataFrame(df)
#print(df)
#f = lambda df:df.loc[df["parameter"].isin(["Hkt_init", 'tropism', 't_mat_P'])]
g = sns.FacetGrid(df.dropna(), col = "sub_readout_name",row ="species",
                    hue = "parameter",
                    sharey = False,
                    sharex = True,
                    height = 2,
                    margin_titles = True,
                    palette="tab10"                    
                    )
g.map(sns.lineplot, "parameter_value","value")
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')

g.add_legend()

#### ScalarMap for colorbar in grey scale ####
color = (0.5,0.5,0.5)

color_start = list(rgb_to_hsv(*color))
color_start[1] = 0
color_start[2] = 1
color_start = hsv_to_rgb(*color_start)

color_end = list(rgb_to_hsv(*color))
color_end[2] = 0.0*color_end[2]
color_end = hsv_to_rgb(*color_end)

cmap = lsc.from_list('ai',[color_start,color,color_end])
m = np.logspace(-.1,.1,10)
norm = plt.Normalize(vmin = np.log10(m[0]),vmax = np.log10(m[-1]))
sm = ScalarMappable(norm=norm, cmap=cmap)

#### 

def create_cm(color):
    color_start = list(rgb_to_hsv(*color))
    color_start[1] = 0.25*color_start[1]
    color_start = hsv_to_rgb(*color_start)

    color_end = list(rgb_to_hsv(*color))
    color_end[2] = .25*color_end[2]
    color_end = hsv_to_rgb(*color_end)

    return lsc.from_list('ai',[color_start,color,color_end])



def my_scatter(*args, data = None,color = None,label = None,):
   
    peaks = data[["parameter_value",'sub_readout_name']]

    df = data.pivot(columns=["readout_name"],index = "parameter_value",values = "value")
    df = df.merge(peaks,on = "parameter_value")
    df.sort_index(inplace=True)

    sns.lineplot(data = df, x = args[0],y =args[1],sort=False, style = "sub_readout_name",
                 color=color)
    
    own_cm = create_cm(color)
    #print(own_cm,color)
    sns.scatterplot(data = df.iloc[::5], x = args[0],y =args[1], style = "sub_readout_name",
                    hue="parameter_value", palette=own_cm,zorder=100, edgecolor='none',
                    s=2,
                    )
#f = lambda df:df.copy()#df.loc[df["species"].isin(["[P]"])]

def my_filter(df, p_type):

    sdf = df.loc[df["readout_name"] == f"{p_type}imum_magnitude"]
    if len(sdf) == 0:
        return pd.DataFrame(columns=df.columns)
    if p_type=='min':
        srn = sdf.loc[sdf["value"].idxmin(),"sub_readout_name"]
    else:
        srn = sdf.loc[sdf["value"].idxmax(),"sub_readout_name"]
    #print(df.loc[df["sub_readout_name"] == srn]["readout_name"].unique())

    return df.loc[df["sub_readout_name"] == srn]

# Function to extract the number after 'peak_'
def extract_peak_number(s):
    match = re.match(r'peak_(\d+)', s)
    if match:
        return int(match.group(1))
    return -1  # Return -1 if not in 'peak_x' format or for handling errors

markers = {
    "peak_0": "o",
    "peak_1": "H",
    "peak_2": ".",
    "peak_3": "^",
    "peak_4": ">",
    "peak_5": "p",
    "peak_6": "o",
    "peak_7": "o",}

# Set rcParams for the desired font sizes
plt.rcParams['axes.titlesize'] = 16      # Title font size
plt.rcParams['axes.labelsize'] = 12      # Axis label font size
plt.rcParams['xtick.labelsize'] = 12     # X tick label font size
plt.rcParams['ytick.labelsize'] = 12     # Y tick label font size
plt.rcParams['legend.fontsize'] = 11     # Legend font size
plt.rcParams['legend.title_fontsize'] = 12 # Legend title font size

for looker in ['min','max']:

    sub_df = df.groupby(df.columns.drop(["readout_name","sub_readout_name","value"],).tolist(),as_index = False).apply(my_filter,looker)

    g = sns.FacetGrid(sub_df.dropna(), col = "species",col_wrap=5,
                        hue = "parameter",
                        sharey = False,
                        sharex = False,
                        height = 1.5,
                        aspect=1,
                        margin_titles = True,
                        palette="tab10",
                        col_order=['[P]','[R]','[E]','parasitemia','[iE]','LDH','RPI','Hb','Hkt','oiE']                    
                        )
    g.map_dataframe(my_scatter,f"{looker}imum_magnitude",f"{looker}imum_timepoint")
    g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
    # Create custom legend
    colors = sns.color_palette("tab10", len(params_bounds.keys()))  # Get 6 colors from the tab10 colormap

    # Find the highest x in 'peak_x' strings using regex
    max_x = max(extract_peak_number(s) for s in df['sub_readout_name'] if re.match(r'peak_\d+', s))

    # Get colors from tab10 colormap
    colors = sns.color_palette("tab10", len(params_bounds.keys()))

    # Create custom legend using labels and markers
    legend_elements = []
    for i, para in enumerate(params_bounds.keys()):
        legend_elements.append(Line2D([0], [0], color=colors[i], label=para))
    i=0
    for k,v in markers.items():
        if i == max_x:break
        label = k
        marker = v  # Default marker is 'o' if not specified in markers dict
        legend_elements.append(Line2D([0], [0], marker=marker, color='grey', label=label,
                                    markerfacecolor=color, markersize=5)
                                )
        i+=1

    # Add the custom legend to the plot
    g.add_legend(title='Parameter', handles=legend_elements, loc='center right', bbox_to_anchor=(1., 0.65))

    cbar_ax = g.figure.add_axes([.9, .2, .01, .22])
    cbar = g.figure.colorbar(sm,cax=cbar_ax, label="Scale of parameter value")

    cbar.locator     = mpl.ticker.LinearLocator(5)
    cbar.formatter   = mpl.ticker.FuncFormatter(lambda x,pos: str(round(10**x,1)))
    cbar.update_ticks()
    cbar.ax.yaxis.set_label_position('right')
    cbar.ax.tick_params(labelsize=11, direction='out', pad=3)

    g.savefig(f"facetgrid_{looker}_timepoint.pdf", dpi=200)

for looker in ['min','max']:

    sub_df = df.groupby(df.columns.drop(["readout_name","sub_readout_name","value"],).tolist(),as_index = False).apply(my_filter,looker)

    g = sns.FacetGrid(sub_df.dropna(), col = "species",col_wrap=5,
                        hue = "parameter",
                        sharey = False,
                        sharex = False,
                        height = 1.7,
                        aspect=1,
                        margin_titles = True,
                        palette="tab10",
                        col_order=['[P]','[R]','[E]','parasitemia','[iE]','LDH','RPI','Hb','Hkt','oiE']                    
                        )
    g.map_dataframe(my_scatter,f"{looker}imum_magnitude",f"{looker}imum_width")
    g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
    # Create custom legend
    colors = sns.color_palette("tab10", len(params_bounds.keys()))  # Get 6 colors from the tab10 colormap

    # Find the highest x in 'peak_x' strings using regex
    max_x = max(extract_peak_number(s) for s in sub_df['sub_readout_name'] if re.match(r'peak_\d+', s))

    # Get colors from tab10 colormap
    colors = sns.color_palette("tab10", len(params_bounds.keys()))

    # Create custom legend using labels and markers
    legend_elements = []
    for i, para in enumerate(params_bounds.keys()):
        legend_elements.append(Line2D([0], [0], color=colors[i], label=para, markerfacecolor=colors[i]))
    
    legend_elements.append(Line2D([0], [0], marker='o', color='w', label='Marker', markerfacecolor='w', markersize=0))
    i=0
    for k,v in markers.items():
        if i == max_x+1:break
        label = k
        marker = v  # Default marker is 'o' if not specified in markers dict
        legend_elements.append(Line2D([0], [0], marker=marker, color='grey', label=label,
                                    markerfacecolor=color, markersize=5)
                                )
        i+=1

    # Add the custom legend to the plot
    g.add_legend(title='Parameter', handles=legend_elements, loc='center right', bbox_to_anchor=(.99, 0.65))

    # not working cbar addition
    cbar_ax = g.figure.add_axes([.9, .2, .01, .22])
    cbar = g.figure.colorbar(sm,cax=cbar_ax, label="Scale of parameter value")

    cbar.locator     = mpl.ticker.LinearLocator(5)
    cbar.formatter   = mpl.ticker.FuncFormatter(lambda x,pos: str(round(10**x,1)))
    cbar.update_ticks()
    cbar.ax.yaxis.set_label_position('right')
    cbar.ax.tick_params(labelsize=11, direction='out', pad=3)

    g.savefig(f"facetgrid_{looker}_width.pdf", dpi=200)

plt.show()
