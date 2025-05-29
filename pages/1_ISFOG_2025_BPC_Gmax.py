import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import os
from insitu import fun_sigeff_v, fun_Gmax_with_depth

st.set_page_config(
    page_title="Sensitivity of your initial Gmax estimate",
    page_icon="📊",
)

st.markdown(
    """
    # Sensitivity of pile monotonic response to your initial Gmax estimate
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

st.write(r"Change the slider for $$G_{ref}$$ and see how the pile's lateral monotonic response changes")


sel_Gref_MPa = st.slider(r'$G_{ref}$', min_value=50., max_value=100., value=65., step=1.)
sel_Gref_kPa = sel_Gref_MPa * 1e3

st.markdown(
    """
    Initial Gmax is calculated as:
    $$G_{max} = G_{ref} \cdot \\frac{(2.17-e_0)^2}{(1+e_0)} \cdot \sqrt{\\frac{p'_0}{100}}$$
    
    $$p'_0 = \\frac{2\cdot K_0 + 1}{3} \cdot \gamma' \cdot z$$
    
    $$e_0 = 0.488$$ and corresponds to $D_R=89 \%$
    
    $$K_0 = 1 - \sin{\phi'_{cs}}$$
    """
)
st.markdown("""---""")
# Add monotonic response here

st.markdown("""---""")
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
Gmax_with_depth_calc = fun_Gmax_with_depth(depths, pit_e0, sel_Gref_MPa, K0, unit_weight)
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
        title_text=r'Initial Gmax (MPa)',
        title_standoff=0,
        range=[0, 150],
        mirror='allticks',
        side='top'
    ),
)

st.markdown("""---""")
st.write(r"Visualise the initial in-situ Gmax predicted by your selection of $G_{ref}$ and compare it to various correlations from the literature")

st.plotly_chart(fig_Gmax, use_container_width=True)
