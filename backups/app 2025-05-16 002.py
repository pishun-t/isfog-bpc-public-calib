import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from insitu import fun_Gmax_with_depth
from elastic import fun_ICG3S_R
from plastic import fun_critical_state_line
from m02integrate import fun_m02_txc
import os

# import plotting data
data_folder = r"data"

# basic properties - constants
pit_depth = 2.6 # metres
depths = np.arange(0, pit_depth, 0.01)
unit_weight = 17.34 # kN/m3
emin, emax = (0.444, 0.788)
pit_e0 = 0.4877
phi_cs = 32.7 # degrees
phi_cs_rad = np.radians(phi_cs)

# streamlit app
st.set_page_config(
    page_title="Calibration - ISFOG 2025 Blind Prediction Contest",
    page_icon=":bar_chart:", # bar chart emoji
    layout="wide",
)

# show sidebar
st.sidebar.title("Choose your own model parameters")

# Add model sliders to sidebar only
st.sidebar.subheader("IC-MAGE-M02 parameters")


st.sidebar.subheader("In-situ condition")


st.title('Calibrate IC-MAGE-M02')
st.write('This app allows you to calibrate model parameters and see the effect on the lateral monotonic pile response.')

st.header("Single element response")

# zeroth section - slides for all material parameters
calib_params = {
    'e_csref': 0.788,
    'l_cs': 0.091,
    'csi_cs': 0.138,
    'pref': 100,
    'mcs': 1.3178,
    'k1': 2.36,
    'k2': 0,
    'l1': 3.42,
    'l2': 0,
    'g0': 65000,
    'mg': 0.5,
    'a0': 5e-5,
    'b': 1.1,
    'rgmin': 0.15,
    'pr': 0.23,
}

# create sliders for selected parameters only
# g0, a0, b, l_cs, csi_cs, k1, l1
col1_ela, col2_pla, col3_txc = st.columns(3)

col1_ela.subheader('Non-linear elastic parameters')
sel_Gref_MPa = col1_ela.slider(r'$G_{ref}$', min_value=50., max_value=100., value=65., step=1., key="Gref_MPa_1")
sel_loga = col1_ela.slider(r'$log_{10}(a)$ (%)', min_value=-3., max_value=0., value=np.log10(5e-3), step=0.01, key="log10_a_1")
sel_b = col1_ela.slider(r'$b$ (slope of stiffness degradation curve)', min_value=0.6, max_value=1.3, value=1.1, step=0.01, key="b_1")

col2_pla.subheader('Plasticity parameters')
sel_l_cs = col2_pla.slider(r'$\lambda$ (Steepness of CSL)', min_value=0.001, max_value=0.2, value=0.091, step=0.001, key="lambda_cs_1")
sel_csi_cs = col2_pla.slider(r"$\xi$ (Curvature of CSL)", min_value=0.001, max_value=0.3, value=0.138, step=0.001, key="xi_cs_1")
sel_k1 = col2_pla.slider(r'$k_1$ (Peak strength)', min_value=0.01, max_value=5., value=2.36, step=0.01, key="k1_1")
sel_l1 = col2_pla.slider(r'$l_1$ (Plastic dilatancy)', min_value=0.01, max_value=5., value=3.42, step=0.01, key="l1_1")

# test inital conditions
col3_txc.subheader('Drained triaxial compression tests')

df_info = pd.read_csv(os.path.join(data_folder, 'lab_mono_init_info.csv'))

df_txc = pd.read_csv(os.path.join(data_folder, 'lab_mono_dtxc_stressstrain_long.csv'))
df_txc = df_txc[~df_txc['test_id'].isin(['CD1'])]

# plot as dashed lines
fig = px.line(df_txc, x='eax_prc', y='q_kPa', color='test_id', line_shape='linear')

fig.update_layout(title="Medium to dense drained triaxial tests",
    yaxis=dict(title_text='q (kPa)',),
    xaxis=dict(title_text=r'Axial strain (%)', range=[0, 20],),
)

# calculate model response based on selected parameters
sel_mat_params = calib_params.copy()

sel_Gref_kPa = sel_Gref_MPa * 1e3
sel_mat_params['g0'] = sel_Gref_kPa
sel_a0 = (10**sel_loga) / 100
sel_mat_params['a0'] = sel_a0
sel_mat_params['b'] = sel_b
sel_mat_params['l_cs'] = sel_l_cs
sel_mat_params['csi_cs'] = sel_csi_cs
sel_mat_params['k1'] = sel_k1
sel_mat_params['l1'] = sel_l1

sel_tests = df_txc['test_id'].unique()

outputs = {}
for test_id in sel_tests:
    info = df_info[df_info['test_id'] == test_id].iloc[0]
    e0, p0, q0 = np.array(info[['e0', 'p0_kPa', 'q0_kPa']], dtype=float)

    output = fun_m02_txc(
        material_parameters=sel_mat_params,
        initial_conditions={'e0': e0, 'p0': p0, 'q0': q0}
    )
    outputs[test_id] = output

    # plot model response
    x_plot = np.array(output['eax']) * 100 # convert to %
    y_plot = np.array(output['q'])
    # find corresponding test_id
    fig.add_trace(
        go.Scatter(
            name=f"Model {test_id}",
            x=x_plot, y=y_plot,
        )
    )

col3_txc.plotly_chart(fig, use_container_width=True)

# streamlit line across the entire page
st.markdown("---")

# Calibration - non linear elasticity section
st.header("Model calibration plots")
st.subheader('Calibration: non-linear elasticity - shear')
st.write('CPT and lab-informed calibration of shear modulus at 1) very small strains, 2) degradation with shear strain')

col1, col2, col3 = st.columns(3)
K0 = 1 - np.sin(phi_cs_rad)

# create sliders for K0, Gref, a and b
sel_K0 = col1.slider(f'$K_0$', min_value=0.3, max_value=0.6, step=0.01)
phi = np.degrees(np.arcsin(1 - sel_K0)) # display equivalent phi'
col1.write(f'(Equivalent $\phi_{{cs}}={phi:.1f}\degree$, assuming $K_0=1-sin(\phi_{{cs}})$)')

# first column
# plot Gmax with depth
# static lines for CPT correlated Gmaxs
df = pd.read_csv(os.path.join(data_folder, 'cpt_corr_Gmax_examples_static_long.csv'))

# sort df by Correlation_method and depth
df = df.sort_values(by=['Correlation_method', 'Depth_m'])

fig_Gmax = px.line(
    df, x="Gmax_MPa", y="Depth_m", color="Correlation_method",
    # line styles dashed
    line_dash="Correlation_method",
    line_dash_sequence=["solid", "dash", "dot", "longdash", "longdashdot"],
)

# recalculate Gmax with depth given slider values
# add static line at Gmax = 65 MPa
Gmax_static = fun_Gmax_with_depth(depths, pit_e0, 65, K0, unit_weight)
fig_Gmax.add_trace(
    go.Scatter(
        x=Gmax_static, y=depths, mode='lines', name=rf'Paper G_max',
        line=dict(color='black', width=4)
    )
)
Gmax_with_depth_calc = fun_Gmax_with_depth(depths, pit_e0, sel_Gref_MPa, sel_K0, unit_weight)
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
static_Gnorm = 65 * fun_ICG3S_R(Ed_plot, 10**5e-3, 1.1, 0.15)
fig_Gnorm.add_trace(
    go.Scatter(
        x=Ed_plot, y=static_Gnorm, mode='lines', name=rf'Paper G_tan',
        line=dict(color='black', width=3)
    )
)
red_factor = fun_ICG3S_R(Ed_plot, 10**sel_loga, sel_b, calib_params['rgmin'])
Gnorm = sel_Gref_MPa * red_factor
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

st.markdown("---")

# CSL

st.subheader('Critical State Line')

col1, col2, col3 = st.columns(3)

# update CSL parameters if so wish?
sel_l_cs_2 = col1.slider(r'$\lambda$ (Steepness of CSL)', min_value=0.001, max_value=0.2, value=0.091, step=0.001, key="lambda_cs_2")

# Synchronize: update shared value if either slider changes
if sel_l_cs != st.session_state["lambda_cs_shared"]:
    st.session_state["lambda_cs_shared"] = sel_l_cs
elif sel_l_cs_2 != st.session_state["lambda_cs_shared"]:
    st.session_state["lambda_cs_shared"] = sel_l_cs_2

# Use the shared value everywhere else
sel_l_cs = st.session_state["lambda_cs_shared"]


p_plot = np.logspace(0, 3, num=100)
# plot CSL - variable
e_plot = fun_critical_state_line(p_plot, calib_params['e_csref'] ,sel_l_cs, sel_csi_cs)
fig_CSL = px.line(
    x=p_plot, y=e_plot,
    labels={'x': 'Mean effective stress (kPa)', 'y': 'Void ratio'},
)
fig_CSL.update_layout(
    title="Critical State Line",
    yaxis=dict(title_text='Void ratio',range=[emin, emax]), # e axis limits
    xaxis=dict(
        title_text='Mean effective stress (kPa)',
        range=[1, 3], # xaxis limits
        type='log'
    ) # use logx
)
col2.plotly_chart(fig_CSL, use_container_width=True)


# # second section
# st.header('Non-linear elasticity - bulk')
# st.write('Default calibration only')

# # fourth section
# st.header('Plasticity - Peak strength')

# # fifth section
# st.header('Plasticity - post-peak dilatancy')

# # sixth section
# st.header('Monotonic single element response')

# # load MathJax
# with open("load-mathjax.js", "r") as f:
#     js = f.read()
#     st.components.v1.html(f"<script>{js}</script>", height=0)