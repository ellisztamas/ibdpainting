import pandas as pd
import numpy as np

import plotly.express as px


def plot_ibd_table(ibd_table:pd.DataFrame, sample_name:str, expected_match:list=[], max_to_plot=10):
    """
    Plot allele sharing across the genome.


    Create a interactive line graph showing genetic distance from a test
    individual to each sample in a panel of reference individuals.

    Parameters
    ==========
    ibd_table: pd.DataFrame
        DataFrame with a row for each window in the genome and a column for each 
        sample in the reference panel. Elements show genetic distance between the 
        test individual and each reference individual in a single window.
        This is generated by ibd_table().
    sample_name: str
        Sample name for the individual to check.
        This must be present in the samples in the input VCF.
    expected_match: list
        List of sample names in the reference panel that are expected to be
        ancestors of the test individual.

    Returns
    =======
    Plotly figure object with subplots for each chromosome, showing window
    position along the x-axis and genetic distance from the test individual to
    each reference sample on the y-axis. Line colour indicates whether a sample
    is an expected parent or not. Rolling over the lines shows which sample is
    which.
    """

    
    # Coerce missing data to NaN for correct column means.
    ibd_table = ibd_table.replace(-9,np.NaN)

    # Identify the candidate names *not* among the top `max_to_plot` columns and remove
    # If `max_to_plot` is less than the number of candidates.
    if max_to_plot < ibd_table.shape[1]-1:
        # Get column-mean IBD for each candidate, allowing for missing data
        ibd_scores_for_each_candidate = np.array(
            [ np.nanmean(ibd_table[col]) for col in ibd_table.keys()[1:] ]
        )
        # Identify the candidate names *not* among the top `max_to_plot` columns
        ix = np.argpartition(ibd_scores_for_each_candidate, max_to_plot)[max_to_plot:] # index positions
        columns_to_drop = ibd_table.keys()[ix+1].to_list() # candidate names to be removed
        # Make sure the expected parents are not among the proscribed
        [ columns_to_drop.remove(x) for x in expected_match if x in columns_to_drop ]        
        ibd_table = ibd_table.drop(columns=columns_to_drop) # drop the candidates

    # Make the table long
    ibd_table = ibd_table.melt(id_vars=['window'], var_name='candidate', value_name='distance')
    # Column indicating which candidates should be plotted a different colour.
    ibd_table['colour'] = np.where(
        ibd_table['candidate'].isin(expected_match), ibd_table['candidate'], "Other"
        )
    # Unique list of labels for the legend, sorted to plot "Other" first
    unique_legend_labels = list(ibd_table['colour'].unique())
    unique_legend_labels.insert(0, unique_legend_labels.pop(unique_legend_labels.index("Other")))

    # Split the 'window' column up into separate columns for chromosome, start and stop positions
    ibd_table[['chr', 'window']] = ibd_table['window'].str.split(":", expand=True)
    ibd_table[['start', 'stop']] = ibd_table['window'].str.split("-", expand=True)
    # start and stop positions should be integers for sensible plotting.
    ibd_table['start'] = ibd_table['start'].astype(int)
    ibd_table['stop'] = ibd_table['stop'].astype(int)
    ibd_table['midpoint'] = (ibd_table['start'] + ibd_table['stop']) / 2


    fig = px.line(
        ibd_table,
        x="midpoint", y="distance", color="colour", line_group="candidate",
        title=sample_name,
        labels={
            'midpoint' : 'Position (bp)',
            'distance' : 'Genetic distance'        
        },
        hover_data=['candidate'],
        color_discrete_sequence=[
                 "gray", "red", "blue"],
        category_orders={'colour': unique_legend_labels},
        facet_row = "chr"
        )
    fig.update_traces(mode="markers+lines")

    return fig


