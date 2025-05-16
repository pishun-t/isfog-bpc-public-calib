import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from lib.insitu import fun_Gmax_with_depth
from lib.elastic import fun_ICG3S_R
import os

# import plotting data
data_folder = r"data"
# constants
pit_depth = 2.6 # metres
depths = np.arange(0, pit_depth, 0.01)
unit_weight = 17.34 # kN/m3
void_ratio = 0.4877
phi_cs = 32.7 # degrees
# original 
K0_orig = 1 - np.sin(np.radians(phi_cs))
Gref_orig = 65 # MPa
loga_orig = np.log10(5e-3) # prc
b_orig = 1.1
Rmin = 0.15

K0 = K0_orig
Gref = Gref_orig
loga = loga_orig
b = b_orig

# streamlit app
st.set_page_config(
    page_title="Calibration - ISFOG 2025 Blind Prediction Contest",
    page_icon=":bar_chart:", # bar chart emoji
    layout="wide",
    # initial_sidebar_state="expanded",
    initial_sidebar_state="auto",
)

st.title('Calibrate IC-MAGE-M02')
st.write('This app allows you to calibrate model parameters and see the effect on the lateral monotonic pile response.')

# first section
st.header('Non-linear elasticity - shear')
st.write('This section allows you to calibrate shear modulus at 1) very small strains, 2) degradation with shear strain')
K0_bounds = [0.3, 0.6]
Gref_bounds = [50, 100]
loga_bounds = [-3., 0.]
b_bounds = [0.6, 1.3]

col1, col2, col3, col4 = st.columns(4)

# create sliders for K0, Gref, a and b
K0 = col1.slider(f'$K_0$', min_value=K0_bounds[0], max_value=K0_bounds[1], value=K0, step=0.01)
phi = np.degrees(np.arcsin(1 - K0)) # display equivalent phi'

col1.write(f'(Equivalent $\phi_{{cs}}={phi:.1f}\degree$, assuming $K_0=1-sin(\phi_{{cs}})$)')
Gref = col2.slider('Gref', min_value=Gref_bounds[0], max_value=Gref_bounds[1], value=Gref, step=1)
col2.write(f"Model $G_{{max}} = G_{{ref}}\cdot f(e) \cdot f(p')$")
col2.write(rf"where $f(e)= \frac{{(2.17-e)^2}}{{(1+e)}}$ and $f(p')= \sqrt{{\frac{{p'}}{{100}}}}$")

loga = col3.slider('log10(a) (%)', min_value=loga_bounds[0], max_value=loga_bounds[1], value=loga, step=0.01)
col3.write(rf"$a={10**loga:.2e}$ (%)")
col3.write(r"$G = G_{max}\cdot R$")
col3.write(r"Reduction factor $R = R_{min} + \frac{1-R_{min}}{1 - ((E_d / a) ^ b)}$")
b = col4.slider('b', min_value=b_bounds[0], max_value=b_bounds[1], value=b, step=0.01)

# create 2 columns
# left = selected Gmax and correlations
# right = shear strain degradation
col1, col2, col3 = st.columns(3)

# first column
# plot Gmax with depth
# static lines for CPT correlated Gmaxs
df = pd.read_csv(os.path.join(data_folder, 'cpt_corr_Gmax_examples_static_long.csv'))

# sort df by Correlation_method and depth
df = df.sort_values(by=['Correlation_method', 'Depth_m'])
# Don't plot the data corresponding to "Rix & Stokoe (1991) - lower bound" and "Rix & Stokoe (1991) - upper bound"
# df = df[~df.Correlation_method.isin(['Rix & Stokoe (1991) - lower bound', 'Rix & Stokoe (1991) - upper bound'])]
fig_Gmax = px.line(
    df, x="Gmax_MPa", y="Depth_m", color="Correlation_method",
    # line styles dashed
    line_dash="Correlation_method",
    line_dash_sequence=["solid", "dash", "dot", "longdash", "longdashdot"],
)

# recalculate Gmax with depth given slider values
# add static line at Gmax = 65 MPa
Gmax_static = fun_Gmax_with_depth(depths, void_ratio, Gref_orig, K0_orig, unit_weight)
fig_Gmax.add_trace(
    go.Scatter(
        x=Gmax_static, y=depths, mode='lines', name=rf'Paper G_max',
        line=dict(color='black', width=4)
    )
)
Gmax_with_depth_calc = fun_Gmax_with_depth(depths, void_ratio, Gref, K0, unit_weight)
fig_Gmax.add_trace(
    go.Scatter(
        x=Gmax_with_depth_calc, y=depths, mode='lines', name=rf'Selected G_max',
        line=dict(color='red', width=4)
    )
)

fig_Gmax.update_layout(
    title="Gmax with depth",
    yaxis=dict(
        title_text='Depth (m)',
        autorange="reversed",
        range=[0, pit_depth], # yaxis limits
    ),
    xaxis=dict(
        title_text=r'Gmax (MPa)',
        title_standoff=0,
        range=[0, 150],
        mirror='allticks',
        side='top'
    ),
)

col2.plotly_chart(fig_Gmax, use_container_width=True)

# second column
df = pd.read_csv(os.path.join(data_folder, 'lab_mono_Gnorm_Ed_long.csv'))
# sort by test_id and Ed_prc
df = df.sort_values(by=['test_id', 'Ed_prc'])

# plot lab tests
fig_Gnorm = px.scatter(
    df, x="Ed_prc", y="Gtan/fe/fp", color="test_id",
)
# ?
Ed_plot = np.logspace(-4,1,num=500)
# plot static line at Gtan = 65 MPa
static_Gnorm = Gref_orig * fun_ICG3S_R(Ed_plot, 10**loga_orig, b_orig, Rmin)
fig_Gnorm.add_trace(
    go.Scatter(
        x=Ed_plot, y=static_Gnorm, mode='lines', name=rf'Paper G_tan',
        line=dict(color='black', width=3)
    )
)
red_factor = fun_ICG3S_R(Ed_plot, 10**loga, b, Rmin)
Gnorm = Gref * red_factor
fig_Gnorm.add_trace(
    go.Scatter(
        x=Ed_plot, y=Gnorm, mode='lines', name=rf'Selected Gtan',
        line=dict(color='red', width=3)
    )
)
fig_Gnorm.update_layout(
    title="Tangent shear modulus Gtan with strain level",
    xaxis=dict(
        title_text=r'Shear strain invariant Ed (%)',
        type='log', # log axis
        range=[-4,0]# range=[1e-4, 1e1], # xaxis limits
    ),
    yaxis=dict(
        title_text = r"G_tan / f(e) / f(p') (MPa)",
        range=[0, 100], # MPa normed
    ),
    legend=dict(title_text='Test ID'), # rename legend title from "test_id" to "Test ID"
)
# add arrow at 86.24 pointing to left at x=1e-4, y=100
fig_Gnorm.add_annotation(
    x=-4, y=86.24, ax=-2, ay=90,
    xref="x", yref="y", axref="x", ayref="y",
    text="Best fit using Rix & Stokoe (1991) mean correlation", align="right",
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black", bgcolor="white",opacity=0.8,
)
fig_Gnorm.add_annotation(
    x=-4, y=85, ax=-2, ay=85,
    xref="x", yref="y", axref="x", ayref="y",
    text="Kokusho (1980) Toyoura sand Gmax", align="right",
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black", bgcolor="white",opacity=0.8,
)
fig_Gnorm.add_annotation(
    x=-4, y=70, ax=-2, ay=70, xref="x", yref="y", axref="x", ayref="y",
    text="Hardin & Richart (1963) Ottawa sand Gmax", align="right",
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black", bgcolor="white",opacity=0.8,
)
fig_Gnorm.add_annotation(
    x=-4, y=65, ax=-2, ay=65, xref="x", yref="y", axref="x", ayref="y",
    text="Winning paper Gmax", align="right",
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black", bgcolor="white",opacity=0.8,
)
col3.plotly_chart(fig_Gnorm, use_container_width=True)

# second section
st.header('Non-linear elasticity - bulk')
st.write('Default calibration only')

# third section
st.header('Critical State Line')

# fourth section
st.header('Plasticity - Peak strength')

# fifth section
st.header('Plasticity - post-peak dilatancy')

# sixth section
st.header('Monotonic single element response')

# load MathJax
with open("load-mathjax.js", "r") as f:
    js = f.read()
    st.components.v1.html(f"<script>{js}</script>", height=0)