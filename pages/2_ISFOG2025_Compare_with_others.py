import streamlit as st
from streamlit_gsheets import GSheetsConnection
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import os
from metamodel import mm_output

conn = st.connection("gsheets", type=GSheetsConnection)
df_saved_submissions = conn.read(usecols=[0, 1, 2, 3])

### --- Data sources --- ###
data_folder = r"data"
df_pit_mono = pd.read_csv(os.path.join(data_folder, 'pit_mono.csv'))

st.write('#### Compare your submission with others! (Monotonic response)')

sim_x = np.array([2., 5., 10., 15., 20., 26., 33., 40., 50., 60.])

plot_dfs = []

# calculate range of responses from participants
for i, row in df_saved_submissions.iterrows():
    Gref_MPa = row['Gref_MPa']
    log10a = row['log10a']
    b = row['b']
    output_Fx = mm_output(sim_x, Gref=Gref_MPa, a=10**log10a, b=b)
    sim_x_plot = [0.] + list(sim_x)  # add zero at the beginning for the plot
    output_Fx = [0.] + list(np.array(output_Fx))  # add zero at the beginning for the plot
    
    plot_df = pd.DataFrame({
        'Name': row['Name'],
        'Ux_mm': sim_x_plot,
        'Fx_kN': output_Fx,
    })
    plot_dfs.append(plot_df)
    
# concatenate all dataframes
df_plot = pd.concat(plot_dfs, ignore_index=True)

fig_mono = px.line(
    df_plot, x="Ux_mm", y="Fx_kN", color="Name",
)

# update all lines to be thin and semi-transparent
fig_mono.update_traces(
    line=dict(width=1),
    opacity=np.min([1.0, 3/len(df_saved_submissions)]),
)

col1, col2, col3 = st.columns(3)
# display your prediction with checkbox
if col1.checkbox('(Spoiler) show current prediction'):
    sel_Gref_MPa = st.session_state['sel_Gref_MPa']  # default value
    sel_loga = st.session_state['sel_loga']  # default value
    sel_b = st.session_state['sel_b']  # default value
    
    model_mono = mm_output(sim_x, Gref=sel_Gref_MPa, a=10**sel_loga, b=sel_b)
    model_mono = [0.] + list(np.array(model_mono))  # add zero at the beginning for the plot
    sim_x_plot = [0.] + list(sim_x)  # add zero at the beginning for the plot
    
    fig_mono.add_trace(
        go.Scatter(
            x=sim_x_plot, y=model_mono,
            mode='lines+markers', name='Your prediction',
            marker=dict(color='red', size=10, symbol='circle'),
        )
    )
    
# display true response with checkbox
if col2.checkbox('Show true response'):
    x_plot = df_pit_mono['Ux_mm']
    y_plot = df_pit_mono['Fx_kN_exp']
    
    fig_mono.add_trace(
        go.Scatter(
            x=x_plot, y=y_plot,
            mode='lines+markers', name='True response',
            marker=dict(color='black', size=10, symbol='square'),
        )
    )

# display winning prediction with checkbox
if col3.checkbox('Show winning prediction'):
    x_plot = df_pit_mono['Ux_mm']
    y_plot = df_pit_mono['Fx_kN_win']
    fig_mono.add_trace(
        go.Scatter(
            x=x_plot, y=y_plot,
            mode='lines+markers', name='Winning blind prediction',
            marker=dict(color='blue', size=3, symbol='diamond'),
        )
    )

# update layout
fig_mono.update_layout(
    # title="",
    yaxis=dict(
        title_text='Lateral force (kN)',
        range=[0, 100],  # y-axis limits
    ),
    xaxis=dict(
        title_text='Lateral displacement (mm)',
        range=[0, 65],
    ),
    legend_title_text='All predictions',
)

# rearrange legend order so that 'Your prediction' comes at the top
fig_mono.update_layout(legend_traceorder="reversed")
    
st.plotly_chart(fig_mono, use_container_width=True)
    
st.write("#### Have a look at everyone's submission in a big table!")
st.dataframe(df_saved_submissions)
