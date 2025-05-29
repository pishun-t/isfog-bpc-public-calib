import streamlit as st
from streamlit_gsheets import GSheetsConnection
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import os
from insitu import fun_Gmax_with_depth
from elastic import fun_ICG3S_R
from metamodel import mm_output

st.set_page_config(
    page_title="Sensitivity of monotonic response to G(Ed)",
    page_icon="📊",
)

st.markdown(
    """
    # Calibrate your Gmax and shear modulus reduction
    """
)

### --- --- ###

### --- Constants --- ###
pit_depth = 2.6 # m
depths = np.arange(0, pit_depth, 0.01)
pit_e0 = 0.4877
unit_weight = 17.34 # kN/m3
K0 = 0.4598

### --- Data sources --- ###
data_folder = r"data"
df_cpt_raw = pd.read_csv(os.path.join(data_folder, 'cpt_raw_data_long.csv'))
df_cpt_gmax = pd.read_csv(os.path.join(data_folder, 'cpt_corr_Gmax_examples_static_long.csv'))
df_Gnorm_lab = pd.read_csv(os.path.join(data_folder, 'lab_mono_Gnorm_Ed_long.csv'))
df_Gnorm_lab = df_Gnorm_lab.sort_values(by=['test_id', 'Ed_prc'])
df_pit_mono = pd.read_csv(os.path.join(data_folder, 'pit_mono.csv'))

### --- Gmax with depth --- ###
st.markdown(
    """
    ##### First, have a look at your initial in-situ Gmax (from $G_{ref}$) and compare it to various correlations from the literature.
    """)

# session state to store the selected Gref value
if 'sel_Gref_MPa' not in st.session_state:
    st.session_state['sel_Gref_MPa'] = 55.0  # default value
# store log10a and b values
if 'sel_loga' not in st.session_state:
    st.session_state['sel_loga'] = -4.9  # default value
if 'sel_b' not in st.session_state:
    st.session_state['sel_b'] = 0.8  # default value

# Slider for Gref
min_max = (55., 95.)
sel_Gref_MPa = st.slider(
    r'$G_{ref}$',
    min_value=min_max[0], max_value=min_max[1],
    step=.1
)
sel_Gref_kPa = sel_Gref_MPa * 1e3
# update session state for Gref
st.session_state['sel_Gref_MPa'] = sel_Gref_MPa

# Explanation
st.markdown(
    """
    Initial Gmax is calculated as:
    $$G_{max} = G_{ref} \cdot \\frac{(2.17-e_0)^2}{(1+e_0)} \cdot \sqrt{\\frac{p'_0}{100}}$$
    
    $$p'_0 = \\frac{2\cdot K_0 + 1}{3} \cdot \gamma' \cdot z$$
    
    $$e_0 = 0.488$$; corresponds to $D_R=89 \%$
    
    $$K_0 = 1 - \sin{\phi'_{cs}} \\approx 0.46$$
    """
)
st.markdown("""---""")

fig_Gmax = px.line(
    df_cpt_gmax, x="Gmax_MPa", y="Depth_m", color="Correlation_method",
    # line styles dashed
    line_dash="Correlation_method",
    line_dash_sequence=["solid", "dash", "dot", "longdash", "longdashdot"],
)

# # recalculate Gmax with depth given slider values
# # add static line at Gmax = 65 MPa
# Gmax_static = fun_Gmax_with_depth(depths, pit_e0, 65, K0, unit_weight)
# fig_Gmax.add_trace(
#     go.Scatter(
#         x=Gmax_static, y=depths, mode='lines', name=rf'Paper G_max',
#         line=dict(color='black', width=4)
#     )
# )
Gmax_with_depth_calc = fun_Gmax_with_depth(depths, pit_e0, sel_Gref_MPa, K0, unit_weight)
fig_Gmax.add_trace(
    go.Scatter(
        x=Gmax_with_depth_calc, y=depths, mode='lines', name='Your selected G_max',
        line=dict(color='red', width=5, dash='dash')
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
        title_text=r'Initial Gmax (MPa)',
        title_standoff=0,
        range=[0, 150],
        mirror='allticks',
        side='top'
    ),
)

st.plotly_chart(fig_Gmax, use_container_width=True)

st.markdown("""---""")

### --- G(Ed) shear strain degradation --- ###
col1, col2 = st.columns(2)
min_max = (-4.9, -3.1)
sel_loga = col1.slider(
    r'$\log_{10}(a)$',
    min_value=min_max[0], max_value=min_max[1],
    step=0.01
)
min_max = (0.8, 1.2)
sel_b = col2.slider(
    r'$b$',
    min_value=min_max[0], max_value=min_max[1],
    step=0.01
)
# update session state for a and b
st.session_state['sel_loga'] = sel_loga
st.session_state['sel_b'] = sel_b

st.markdown(
    """
    ##### Then, try to fit the normalised shear modulus degradation curve against the available laboratory data at medium-to-large shear strains.
    """
)
st.markdown(
    """
    $$G_{tan}= G_{max} \cdot (R_G + \\frac{1-R_G}{1 + (E_d / a) ^b})$$
    
    $R_G = 0.15$ is the minimum reduction factor for the elastic tangent shear modulus;
    $a, b$ are your selected parameters from the sliders above.        
    """
)

fig_Gnorm = px.scatter(
    df_Gnorm_lab, x="Ed_prc", y="Gtan/fe/fp", color="test_id",
)

Ed_plot = np.logspace(-4,-1,num=200)

# add horizontal line at sel_Gref_MPa
fig_Gnorm.add_trace(
    go.Scatter(
        x=Ed_plot, y=np.full_like(Ed_plot, sel_Gref_MPa), mode='lines', name='Your Gref',
        line=dict(color='black', width=2, dash='dashdot')
    )
)

# # plot static line at Gtan = 65 MPa
# static_Gnorm = 65 * fun_ICG3S_R(Ed_plot, calib_params['a0']*100, calib_params['b'], calib_params['rgmin'])
# fig_Gnorm.add_trace(
#     go.Scatter(
#         x=Ed_plot, y=static_Gnorm, mode='lines', name=rf'Paper Gtan',
#         line=dict(color='black', width=3)
#     )
# )
red_factor = fun_ICG3S_R(Ed_plot, (10**sel_loga) * 100, sel_b, 0.15)
Gnorm = sel_Gref_MPa * red_factor
fig_Gnorm.add_trace(
    go.Scatter(
        x=Ed_plot, y=Gnorm, mode='lines', name=r'Your best a, b',
        line=dict(color='red', width=5, dash='dash')
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

# # add arrow at 86.24 pointing to left at x=1e-4, y=100
# fig_Gnorm.add_annotation(
#     x=-4, y=86.24, ax=-2, ay=90,
#     xref="x", yref="y", axref="x", ayref="y",
#     text="Best fit using Rix & Stokoe (1991) mean correlation", align="right",
#     showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
#     # bgcolor="white",opacity=0.8,
# )
# fig_Gnorm.add_annotation(
#     x=-4, y=85, ax=-2, ay=85,
#     xref="x", yref="y", axref="x", ayref="y",
#     text="Kokusho (1980) Toyoura sand Gmax", align="right",
#     showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
#     # bgcolor="white",opacity=0.8,
# )
# fig_Gnorm.add_annotation(
#     x=-4, y=70, ax=-2, ay=70, xref="x", yref="y", axref="x", ayref="y",
#     text="Hardin & Richart (1963) Ottawa sand Gmax", align="right",
#     showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
#     # bgcolor="white",opacity=0.8,
# )
# fig_Gnorm.add_annotation(
#     x=-4, y=65, ax=-2, ay=65, xref="x", yref="y", axref="x", ayref="y",
#     text="Winning paper Gmax", align="right",
#     showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black",
#     # bgcolor="white",opacity=0.8,
# )
st.plotly_chart(fig_Gnorm, use_container_width=True)

st.markdown("""---""")

st.markdown("""
    ##### Submit and go to the next page to compare your results with others!
""")

# submit button
st.write(f"Your best parameters are: $G_{{ref}}={sel_Gref_MPa:.1f}$ MPa, $log10a={sel_loga:.2f}$, $b={sel_b:.2f}$")

# Add name input field
submission_name = st.text_input("Please enter your unique name (doesn't have to be real):", "")


conn = st.connection("gsheets", type=GSheetsConnection)
df = conn.read(worksheet="ISFOG2025_Compare_with_others", usecols=[0, 1, 2, 3])

df_new = pd.concat([
    df,
    pd.DataFrame({
        'Name': [submission_name],
        'Gref_MPa': [sel_Gref_MPa],
        'log10a': [sel_loga],
        'b': [sel_b],
    })
    ])

if st.button("Submit your parameters!"):
    
    # check all fields have been filled
    if submission_name == "":
        st.error("Please enter your name!")
    elif submission_name in list(df['Name']):
        st.error("This name has already been used! Please choose a different name.")
    elif sel_Gref_MPa is None or sel_loga is None or sel_b is None:
        st.error("Please select all parameters!")
    else:
        df_test = conn.update(
            data=df_new
        )
        st.cache_data.clear()
        st.success(f'Hi {submission_name}, your parameters have been submitted!')

# Plot
sim_x = np.array([2., 5., 10., 15., 20., 26., 33., 40., 50., 60.])

model_mono = mm_output(sim_x, Gref=sel_Gref_MPa, a=10**sel_loga, b=sel_b)
model_mono = [0.] + list(np.array(model_mono))
sim_x = [0.] + list(sim_x)  # add zero at the beginning for the plot

fig_mono_response = go.Figure()
fig_mono_response.add_trace(
    go.Scatter(
        x=sim_x, y=model_mono,
        mode='lines+markers', name='Your prediction!',
        # no line, just markers
        marker=dict(color='red', size=10, symbol='circle'),
    )
)

fig_mono_response.update_layout(
    # show legend
    showlegend=True,
    title=f"Predicted monotonic response from {submission_name}",
    yaxis=dict(
        title_text='Lateral force (kN)',
        range=[0, 100], # yaxis limits
    ),
    xaxis=dict(
        title_text='Lateral displacement (mm)',
        range=[0, 65],
    ),
)

st.plotly_chart(fig_mono_response, use_container_width=True)
