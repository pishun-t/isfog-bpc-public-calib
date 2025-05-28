import streamlit as st

st.set_page_config(
    page_title="M02 calibration",
    page_icon=":bar_chart:", # bar chart emoji
)

st.write("# Welcome to the IC-MAGE-M02 calibration web-app")

st.markdown(
    """
    This web-app uses the TU Darmstadt blind prediction contest winning method to
    illustrate the calibration of the IC-MAGE-M02 model parameters.
    The model is non-linear elastic with a Mohr-Coulomb based plastic yield surface dependent on 
    Been & Jefferies (1985) sand state parameter.
    The app can take you through the calibration of parameters from
     in-situ information, laboratory data, the single element response, all the way to 
      the eventual effect on the lateral monotonic pile response.
"""
)








import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from insitu import fun_Gmax_with_depth
from elastic import fun_ICG3S_R
from plastic import fun_critical_state_line
from m02integrate import fun_m02_txc
import os












##### ----- DATA SOURCES ----- #####
data_folder = r"data"
df_cpt_raw = pd.read_csv(os.path.join(data_folder, 'cpt_raw_data_long.csv'))
df_info = pd.read_csv(os.path.join(data_folder, 'lab_mono_init_info.csv'))
df_txc = pd.read_csv(os.path.join(data_folder, 'lab_mono_dtxc_stressstrain_long.csv'))
df_txc = df_txc[~df_txc['test_id'].isin(['CD1'])]
df_cpt_gmax = pd.read_csv(os.path.join(data_folder, 'cpt_corr_Gmax_examples_static_long.csv'))
df_cpt_gmax = df_cpt_gmax.sort_values(by=['Correlation_method', 'Depth_m'])
df_Gnorm_lab = pd.read_csv(os.path.join(data_folder, 'lab_mono_Gnorm_Ed_long.csv'))
df_Gnorm_lab = df_Gnorm_lab.sort_values(by=['test_id', 'Ed_prc'])
df_csl_end = pd.read_csv(os.path.join(data_folder, 'lab_csl_endpoints.csv'))
df_csl_ep_dtxc = pd.read_csv(os.path.join(data_folder, 'lab_csl_dtxc_ep_long.csv'))
df_plastic = pd.read_csv(os.path.join(data_folder, 'lab_mono_dtxc_plast_clean.csv'))

##### ----- Basic properties - constants ----- #####
pit_depth = 2.6 # metres
depths = np.arange(0, pit_depth, 0.01)
unit_weight = 17.34 # kN/m3
emin, emax = (0.444, 0.788)
pit_e0 = 0.4877
phi_cs = 32.7 # degrees
phi_cs_rad = np.radians(phi_cs)

##### ----- Paper calibration parameters ----- #####
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

##### Initialise streamlit app #####
st.set_page_config(
    page_title="Calibration - ISFOG 2025 Blind Prediction Contest",
    page_icon=":bar_chart:", # bar chart emoji
    layout="wide",
)

# show sidebar
sb = st.sidebar
sb.title("Choose your own IC-MAGE-M02 parameters")

# Add model sliders to sidebar only
sb.subheader("Non-linear elastic parameters")

sel_Gref_MPa = sb.slider(r'$G_{ref}$', min_value=50., max_value=100., value=65., step=1.)
sel_loga = sb.slider(r'$log_{10}(a)$ (%)', min_value=-3., max_value=0., value=np.log10(5e-3), step=0.01)
sel_b = sb.slider(r'$b$ (slope of stiffness degradation curve)', min_value=0.6, max_value=1.3, value=1.1, step=0.01)

sb.subheader('Plasticity parameters')
sel_lambda_cs = sb.slider(r'$\lambda$ (Steepness of CSL)', min_value=0.001, max_value=0.2, value=0.091, step=0.001)
sel_csi_cs = sb.slider(r"$\xi$ (Curvature of CSL)", min_value=0.001, max_value=0.3, value=0.138, step=0.001)
sel_k1 = sb.slider(r'$k_1$ (Peak strength)', min_value=0.01, max_value=5., value=2.36, step=0.01)
sel_l1 = sb.slider(r'$l_1$ (Plastic dilatancy)', min_value=0.01, max_value=5., value=3.42, step=0.01)

sb.subheader("In-situ condition")
sel_K0 = sb.slider(f'$K_0$', min_value=0.3, max_value=0.6, step=0.01)
phi = np.degrees(np.arcsin(1 - sel_K0)) # display equivalent phi'
# create sliders for K0, Gref, a and b
sb.write(f'(Equivalent to $\phi_{{cs}}={phi:.1f}\degree$, assuming $K_0=1-sin(\phi_{{cs}})$)')
K0 = 1 - np.sin(phi_cs_rad)

# Unit conversions
sel_Gref_kPa = sel_Gref_MPa * 1e3
sel_a0 = (10**sel_loga) / 100

# model specific
sel_mat_params = calib_params.copy()
sel_mat_params['g0'] = sel_Gref_kPa
sel_mat_params['a0'] = sel_a0
sel_mat_params['b'] = sel_b
sel_mat_params['l_cs'] = sel_lambda_cs
sel_mat_params['csi_cs'] = sel_csi_cs
sel_mat_params['k1'] = sel_k1
sel_mat_params['l1'] = sel_l1

st.title('Calibrate IC-MAGE-M02')
st.write('This app allows you to calibrate model parameters and see the effect on the lateral monotonic pile response.')

# PART 1
st.header("Single element response")

col1_txc, col2_txc = st.columns(2)

col1_txc.subheader('Medium to dense drained triaxial compression tests')
col2_txc.subheader('') # spacer only

# plot lab tests - stress-strain
fig1 = px.line(df_txc, x='eax_prc', y='q_kPa', color='test_id', line_shape='linear')
fig1.update_layout(title="Stress-strain",
    yaxis=dict(title_text='q (kPa)',),
    xaxis=dict(title_text=r'Axial strain (%)', range=[0, 20],),
)
fig2 = px.line(df_txc, x='eax_prc', y='evol_prc', color='test_id', line_shape='linear')
fig2.update_layout(title="Shear-volumetric coupling",
    yaxis=dict(title_text=r'Vol. strain (%)',),
    xaxis=dict(title_text=r'Axial strain (%)', range=[0, 20],),
)

# Model response
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
    # find corresponding test_id
    fig1.add_trace(
        go.Scatter(name=f"Model", x=output['eax_prc'], y=output['q'],
                   line=dict(color='black', width=2, dash='dash'))
    )
    fig2.add_trace(
        go.Scatter(name=f"Model", x=output['eax_prc'], y=output['evol_prc'],
                   line=dict(color='black', width=2, dash='dashdot'))
    )

col1_txc.plotly_chart(fig1, use_container_width=True)
col2_txc.plotly_chart(fig2, use_container_width=True)

# streamlit line across the entire page
st.markdown("---")

# Calibration - non linear elasticity section
st.header("Model calibration plots")
st.subheader('Calibration: non-linear elasticity (shear)')
st.write('CPT and lab-informed calibration of shear modulus at 1) very small strains, 2) degradation with shear strain')

col1, col2 = st.columns(2)

# first column - plot qc values with depth - all raw data
fig_qc = px.line(df_cpt_raw[df_cpt_raw.value_type=='qc'], x='value', y='Depth_m', color='data_location_id')
# for the "ave" label, change colour to black
fig_qc.update_traces(selector=dict(name='ave'), line_color='black', line_width=3)
fig_qc.update_layout(
    title='Cone resistance',
    yaxis=dict(
        title_text='Depth (m)',
        autorange="reversed",
        range=[0, pit_depth], # yaxis limits
    ),
    xaxis=dict(
        title_text='qc (MPa)',
        title_standoff=0,
        mirror='allticks',
        side='top',
        range=[0,50],
    ),
    legend=dict(title_text='Location'), # rename legend title from "data_location_id" to "Location"
)

col1.plotly_chart(fig_qc, use_container_width=True)

# plot Gmax with depth
fig_Gmax = px.line(
    df_cpt_gmax, x="Gmax_MPa", y="Depth_m", color="Correlation_method",
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
    title="In-situ Gmax",
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

st.markdown("---")

col1, col2 = st.columns(2)
# plot lab tests
fig_Gnorm = px.scatter(
    df_Gnorm_lab, x="Ed_prc", y="Gtan/fe/fp", color="test_id",
)

Ed_plot = np.logspace(-4,1,num=500)

# plot static line at Gtan = 65 MPa
static_Gnorm = 65 * fun_ICG3S_R(Ed_plot, calib_params['a0']*100, calib_params['b'], calib_params['rgmin'])
fig_Gnorm.add_trace(
    go.Scatter(
        x=Ed_plot, y=static_Gnorm, mode='lines', name=rf'Paper Gtan',
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
    title="Fit Gtan degradation curve",
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
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
    # bgcolor="white",opacity=0.8,
)
fig_Gnorm.add_annotation(
    x=-4, y=85, ax=-2, ay=85,
    xref="x", yref="y", axref="x", ayref="y",
    text="Kokusho (1980) Toyoura sand Gmax", align="right",
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
    # bgcolor="white",opacity=0.8,
)
fig_Gnorm.add_annotation(
    x=-4, y=70, ax=-2, ay=70, xref="x", yref="y", axref="x", ayref="y",
    text="Hardin & Richart (1963) Ottawa sand Gmax", align="right",
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
    # bgcolor="white",opacity=0.8,
)
fig_Gnorm.add_annotation(
    x=-4, y=65, ax=-2, ay=65, xref="x", yref="y", axref="x", ayref="y",
    text="Winning paper Gmax", align="right",
    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
    # bgcolor="white",opacity=0.8,
)
col1.plotly_chart(fig_Gnorm, use_container_width=True)

# plot model response
fig_modelG = go.Figure()

for test_id in sel_tests:
    output = outputs[test_id]
    # plot model response
    fig_modelG.add_trace(
        go.Scatter(name=f"{test_id}", x=output['ed_prc'] * np.sqrt(3), y=output['gtan']/1000,
                   # don't plot line, just points
                    mode='markers', marker=dict(size=5, symbol='circle'),
        )
                #    line=dict(color='black', width=2, dash='dash')) # no line, just points
                #    
    )

fig_modelG.update_layout(
    title="Single element model response",
    xaxis=dict(
        title_text=r'Shear strain invariant Ed (%)',
        type='log', # log axis
        range=[-4,0]# range=[1e-4, 1e1], # xaxis limits
    ),
    yaxis=dict(
        title_text = r"Gtan (MPa)",
        range=[0, 200], # MPa normed
    ),
    legend=dict(title_text='Simulated test ID'), # rename legend title from "test_id" to "Test ID"
)
col2.plotly_chart(fig_modelG, use_container_width=True)

st.markdown("---")

##### ---- Critical state line section ----- #####
st.subheader('Critical State Line')

col1, col2 = st.columns(2)

p_plot = np.logspace(0, 3, num=100)

# Plot static line
# fig_CSL = go.Figure()
fig_CSL = px.line(
    df_csl_ep_dtxc, x="p_kPa", y="void_ratio", color="test_id",
)

subdf_ciu_true_end = df_csl_end[df_csl_end.test_id.str.contains('CU')] # & df_csl_end.type == 'actual']
# print(subdf_cid_true_end)
fig_CSL.add_trace(go.Scatter(
        name="CIU test end points",
        x=subdf_ciu_true_end['p_kPa_end'],
        y=subdf_ciu_true_end['e_end'],
        mode='markers',
        marker=dict(size=10, symbol='cross', color='purple', opacity=0.5),
))
# add horizontal line at e=pit_e0
fig_CSL.add_hline(y=pit_e0, line_color='black', line_width=2, line_dash='dash',
              annotation_text="In-situ e0", 
              annotation_position="bottom left")

fig_CSL.add_trace(
    go.Scatter(
        name=f"Paper CSL\n lambda={calib_params['l_cs']}, xi={calib_params['csi_cs']}",
        x=p_plot,
        y=fun_critical_state_line(p_plot, calib_params['e_csref'], calib_params['l_cs'], calib_params['csi_cs']),
        # color='black',
        line=dict(color='black', width=3),
    )
)

# update line colour to black
# fig_CSL.update_traces(line_color='black', line_width=3)
e_plot = fun_critical_state_line(p_plot, calib_params['e_csref'] ,sel_lambda_cs, sel_csi_cs)

# Plot selected line
fig_CSL.add_trace(
    go.Scatter(name=f"Selected CSL", x=p_plot, y=e_plot,
               line=dict(color='red', width=3))
)

# Plot end points

fig_CSL.update_layout(
    title="Critical State Line",
    yaxis=dict(title_text='Void ratio',range=[emin, emax]), # e axis limits
    xaxis=dict(
        title_text='Mean effective stress (kPa)',
        range=[1, 3], # xaxis limits
        type='log'
    ) # use logx
)
col1.plotly_chart(fig_CSL, use_container_width=True)

st.markdown("---")
col1, col2 = st.columns(2)
col1.header('Plasticity - peak strength')

e_peaks = df_plastic['e_peak'].to_numpy()
Mcs = calib_params['mcs']
Mpeak = df_plastic['Mpeak'].to_numpy()
Mpeaks_Mcs = Mpeak - Mcs
p_peaks_assumed =  3 * df_info.loc[df_info.test_id.str.contains('CD'), "p0_kPa"].to_numpy() / (3-Mpeak)

# plot peak strength
psi_peaks = e_peaks - fun_critical_state_line(p_peaks_assumed, calib_params['e_csref'], sel_lambda_cs, sel_csi_cs)

fig_peak = px.scatter(
    x=psi_peaks, y=Mpeaks_Mcs, color=df_plastic['test_id'],
    title="Peak strength",
    labels={'x': 'psi State parameter', 'y': 'M - Mcs'},
)
x_plot_psi = np.linspace(np.min(psi_peaks), 0.0, num=2)
y_plot = -sel_k1 * x_plot_psi
fig_peak.add_trace(
    go.Scatter(name=f"Selected peak strength",
        x=x_plot_psi,
        y=y_plot,
        line=dict(color='red', width=3),
    )
)
col1.plotly_chart(fig_peak, use_container_width=True)

# plot plastic dilatancy
col2.header('Plasticity - post-peak dilatancy')