import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from lib.insitu import fun_Gmax_with_depth
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

K0 = 1 - np.sin(np.radians(phi_cs))
Gref = 65

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

col1, col2 = st.columns(2)
# create sliders for K0 and Gref
K0 = col1.slider(f'$K_0$', min_value=K0_bounds[0], max_value=K0_bounds[1], value=K0, step=0.01)
phi = np.degrees(np.arcsin(1 - K0)) # display equivalent phi'

col1.write(f'(Equivalent $\phi_{{cs}}={phi:.1f}\degree$, assuming $K_0=1-sin(\phi_{{cs}})$)')
Gref = col2.slider('Gref', min_value=Gref_bounds[0], max_value=Gref_bounds[1], value=Gref, step=1)
col2.write(f"Model $G_{{max}} = G_{{ref}}\cdot f(e) \cdot f(p')$")
col2.write(rf"where $f(e)= \frac{{(2.17-e)^2}}{{(1+e)}}$ and $f(p')= \sqrt{{\frac{{p'}}{{100}}}}$")

# add a horizontal line across the screen
st.markdown(
    """
    <style>
    .horizontal-line {
        border: 1px solid #000;
        margin: 20px 0;
    }
    </style>
    <div class="horizontal-line"></div>
    """,
    unsafe_allow_html=True
)
# create 4 columns
col1, col2, col3, col4 = st.columns(4)

# first column - plot qc values with depth - all raw data
df = pd.read_csv(os.path.join(data_folder, 'cpt_raw_data_long.csv'))

# plot 
fig_qc = px.line(df[df.value_type=='qc'], x='value', y='Depth_m', color='data_location_id')
# for the "ave" label, change colour to black
fig_qc.update_traces(selector=dict(name='ave'), line_color='black', line_width=3)
fig_qc.update_layout(
    title='qc with depth',
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
    ),
    legend=dict(title_text='Location'), # rename legend title from "data_location_id" to "Location"
)

col1.plotly_chart(fig_qc, use_container_width=True)

# second column - plot fs values with depth - all raw data
fig_fs = px.line(df[df.value_type=='fs'], x='value', y='Depth_m', color='data_location_id')
# rename legend title from "data_location_id" to "Location"
fig_fs.for_each_trace(
    lambda t: t.update(
        name=t.name.split('=')[1] if '=' in t.name else t.name,
        legendgroup=t.name.split('=')[1] if '=' in t.name else t.name
    )
)
# for the "ave" label, change colour to black
fig_fs.update_traces(selector=dict(name='ave'), line_color='black', line_width=3)
fig_fs.update_layout(
    title="fs with depth",
    yaxis=dict(
        title_text='Depth (m)',
        autorange="reversed",
        range=[0, pit_depth], # yaxis limits
    ),
    xaxis=dict(
        title_text='fs (MPa)',
        title_standoff=0,
        mirror='allticks',
        side='top'
    ),
    legend=dict(title_text='Location'), # rename legend title from "data_location_id" to "Location"
)
col2.plotly_chart(fig_fs, use_container_width=True)

# third column
# plot Gmax with depth
# fig_Gmax = go.Figure()
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
        title_text='fs (MPa)',
        title_standoff=0,
        range=[0, 150],
        mirror='allticks',
        side='top'
    ),
)

col3.plotly_chart(fig_Gmax, use_container_width=True)



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

