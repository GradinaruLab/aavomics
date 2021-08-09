import numpy
import pandas
import anndata

from scipy import signal
from enum import Enum
from diffxpy import api as diffxpy
from scipy import stats

from plotly import graph_objects
from plotly import offline as plotly


def get_background_trough(total_transcript_counts):

    bin_range = list(range(1, max(total_transcript_counts.astype(numpy.int32))+1))

    counts, bin_edges = numpy.histogram(total_transcript_counts, bins=bin_range)

    peaks, prominences = signal.find_peaks(counts, prominence=0)

    first_peak = bin_edges[peaks[numpy.argmax(prominences["prominences"])]]

    negative_counts = -counts[0:first_peak]

    trough_range = list(range(1, first_peak))

    for width in trough_range:

        peaks, _ = signal.find_peaks(negative_counts, width=width)

        if len(peaks) == 1:
            break

    try:
        trough = bin_edges[peaks[0]]
    except:

        peaks, prominences = signal.find_peaks(negative_counts, prominence=0)

        trough = bin_edges[peaks[numpy.argmax(prominences["prominences"])]]
        
    return trough

def plot_clusters(
    coordinates,
    clusters,
    filename=None,
    text=None
):
    coordinates = numpy.array(coordinates)
    clusters = numpy.array(clusters)
    
    traces = []

    for cluster in numpy.unique(clusters):
        
        scatter_trace = graph_objects.Scattergl(
            x=coordinates[clusters == cluster, 0],
            y=coordinates[clusters == cluster, 1],
            mode="markers",
            marker={
                "size": 3
            },
            name=str(cluster),
            text=text
        )
        
        traces.append(scatter_trace)
    
    layout = {}
    
#     layout["height"] = 800
#     layout["width"] = 800
    layout["plot_bgcolor"] = "rgba(255, 255, 255, 0)"
    layout["paper_bgcolor"] = "rgba(255, 255, 255, 0)"
    
    layout = graph_objects.Layout(layout)
    
    figure = graph_objects.Figure(
        data=traces,
        layout=layout
    )
    
    if filename is not None:
        if filename.endswith(".html"):
            figure.write_html(filename)
        else:
            figure.write_image(filename)
    else:
        plotly.iplot(figure)

def plot_gene_expression(
    coordinates,
    gene_counts,
    filename=None,
    text=None
):
    
    coordinates = numpy.array(coordinates)
    gene_counts = numpy.array(gene_counts).reshape((-1, ))
    
    if text is None:
        text = gene_counts
    else:
        text = ["%s %.2f" % (gene_count, text) for gene_count, text in zip(gene_counts, text)]
    
    traces = []

    scatter_trace = graph_objects.Scattergl(
        x=coordinates[:, 0],
        y=coordinates[:, 1],
        marker={
            "color": gene_counts,
            "showscale": True,
            "size": 3
        },
        mode="markers",
        text=text
    )

    traces.append(scatter_trace)
    
    layout = {}
    
#     layout["height"] = 800
#     layout["width"] = 800
    layout["plot_bgcolor"] = "rgba(255, 255, 255, 0)"
    layout["paper_bgcolor"] = "rgba(255, 255, 255, 0)"
    
    layout = graph_objects.Layout(layout)
    
    figure = graph_objects.Figure(
        data=traces,
        layout=layout
    )
    
    if filename is not None:
        figure.write_html(filename)
    else:
        plotly.iplot(figure)


def plot_histogram(values,
                   num_bins=None,
                  is_integer_values=None):

    if is_integer_values is None:
        is_integer_values = True

        # There has got to be a better way to do this...
        # Sample 10% of the values to see if they're all integers
        num_values_to_sample = int(numpy.ceil(numpy.sqrt(len(values))))
        for i in range(num_values_to_sample):
            random_index = numpy.random.randint(0, len(values))
            if values[random_index] != int(values[random_index]):
                is_integer_values = False
                break

    figure_traces = []

    if num_bins is None:
        # Estimate the number of bins using the Rice Rule
        num_bins = int(2 * numpy.ceil(numpy.power(len(values), 1 / 3)))

        if is_integer_values and num_bins > max(values):
            num_bins = int(numpy.ceil(max(values)))

    if is_integer_values:

        min_value = min(values)
        max_value = max(values)

        value_range = max_value - min_value
        bin_size = numpy.ceil(value_range / num_bins)
        min_bin = numpy.floor(min_value / bin_size) * bin_size
        max_bin = numpy.ceil((max_value+1) / bin_size) * bin_size
        bins = numpy.arange(min_bin, max_bin+bin_size, bin_size)
    else:
        bins=num_bins

    y, x = numpy.histogram(values, bins=bins)

    histogram = graph_objects.Bar(
        x=x,
        y=y,
    )

    figure_traces.append(histogram)

    layout_parameters = {
        "hovermode": "closest",
        "bargap": 0
    }
    
    layout_parameters["height"] = 800
    layout_parameters["width"] = 800
    layout_parameters["plot_bgcolor"] = "rgba(255, 255, 255, 0)"
    layout_parameters["paper_bgcolor"] = "rgba(255, 255, 255, 0)"
    
    layout = graph_objects.Layout(layout_parameters)

    figure = graph_objects.Figure(data=figure_traces, layout=layout)

    plotly.iplot(figure)
    

def get_neg_binomial_diffxpy(values):

    X = numpy.array([values for _ in range(2)]).transpose()
    
    obs_df = pandas.DataFrame(index=[str(x) for x in range(X.shape[0])])
    var_df = pandas.DataFrame(index=["0", "1"])

    adata = anndata.AnnData(X=X, obs=obs_df, var=var_df)

    x = diffxpy.fit.model(data=adata, formula_loc="~ 1", batch_size=(int(1e9), 1))

    a = x.a_var[0][0]
    b = x.b_var[0][0]

    n = numpy.exp(b)
    p = 1 - numpy.exp(a)/(numpy.exp(a) + numpy.exp(b))
    
    return scipy.stats.nbinom(n, p)
    

class Infection_Rate_Method(Enum):
    THRESHOLD = 1
    KS = 2
    COUNTING = 3
    CUMULATIVE_PROBABILITY = 4
    RELATIVE_ECDF_AREA = 5
    PAIRWISE_COMPARISONS = 6
    NBINOM_P_NONZERO = 7
    NBINOM_FIXED_P = 8


def shoelace_formula(points):
    
    area = 0
    
    for vertex_index in range(len(points) - 1):
        area += points[vertex_index][0] * points[vertex_index + 1][1]
        area -= points[vertex_index + 1][0] * points[vertex_index][1]
    
    area += points[-1][0] * points[0][1]
    area -= points[0][0] * points[-1][1]
    
    return 0.5 * numpy.abs(area)


def get_transcript_presence_rate(
    signal_transcript_counts,
    method=Infection_Rate_Method.THRESHOLD,
    background_transcript_counts = None,
    threshold=None,
    resolution=100,
    round_to_num_places=None
):
        
    if round_to_num_places is not None:
        signal_transcript_counts = numpy.round(signal_transcript_counts, round_to_num_places)
        background_transcript_counts = numpy.round(background_transcript_counts, round_to_num_places)
    
    if method == Infection_Rate_Method.THRESHOLD:
        
        if threshold is None:
            raise ValueError("Must specify threshold for threshold method")
        
        return (signal_transcript_counts >= threshold).sum()/signal_transcript_counts.shape[0]
    
    elif method == Infection_Rate_Method.KS:
        
        d, _ = scipy.stats.ks_2samp(background_transcript_counts, signal_transcript_counts)
        
        return d
    
    elif method == Infection_Rate_Method.COUNTING:
        
        p_values = numpy.linspace(0, 1, resolution+1)

        unique_noise_counts = set(background_transcript_counts)
        unique_signal_counts = set(signal_transcript_counts)

        unique_counts = sorted(list(unique_noise_counts.union(unique_signal_counts)))

        signal_count_counts = {}
        noise_count_counts = {}

        for count in unique_counts:
            signal_count_counts[count] = signal_transcript_counts[signal_transcript_counts == count].shape[0]
            noise_count_counts[count] = background_transcript_counts[background_transcript_counts == count].shape[0]

        scores = []

        for p_debris in p_values:

            signal_unaccounted_for = 0
            noise_unaccounted_for = 0
            total_noise_estimate = 0
            total_signal_estimate = 0

            for count in unique_counts:

                noise_count = (noise_count_counts[count] * p_debris) /\
                    background_transcript_counts.shape[0]*signal_transcript_counts.shape[0] 
                signal_count = signal_count_counts[count]

                # Log any unaccounted for noise, if it exists
                noise_unaccounted_for += max(0, noise_count - (signal_count + signal_unaccounted_for))

                # Log any unaccounted for signal, if it exists
                signal_unaccounted_for = max(0, (signal_count + signal_unaccounted_for) - noise_count)

                total_noise_estimate += noise_count
                total_signal_estimate += signal_count

            if total_noise_estimate == 0:
                score = signal_unaccounted_for/total_signal_estimate/2
            else:
                score = (signal_unaccounted_for/total_signal_estimate+noise_unaccounted_for/total_noise_estimate)/2
                
            scores.append(score)
          
        estimated_infectivity_rate = 1-p_values[numpy.argmin(scores)]
        
        return estimated_infectivity_rate, scores
    
    elif method == Infection_Rate_Method.CUMULATIVE_PROBABILITY:
        
        p_values = numpy.linspace(0, 1, resolution+1)
    
        max_count = max(background_transcript_counts)

        unique_counts = set(background_transcript_counts).union(set(signal_transcript_counts)).difference([0])
        unique_counts = sorted(list(unique_counts))

        distance_values = []

        for p_background in p_values:

            distance = 0

            previous_count = 0
            previous_p_x_1 = background_transcript_counts[background_transcript_counts == 0].shape[0] * \
                p_background/background_transcript_counts.shape[0]
            previous_p_x_2 = signal_transcript_counts[signal_transcript_counts == 0].shape[0] / \
                signal_transcript_counts.shape[0]

            for count in unique_counts:

                if count > max_count:
                    break

                # For this transcript count, calculate the expected probability of seeing it, assuming our debris percentage
                p_x_1 = background_transcript_counts[background_transcript_counts <= count].shape[0] * \
                    p_background/background_transcript_counts.shape[0]

                # Similarly, calculate the probability of this transcript count in the debris + cell data
                p_x_2 = signal_transcript_counts[signal_transcript_counts <= count].shape[0]/signal_transcript_counts.shape[0]

                distance += numpy.abs(p_x_1-p_x_2)

            distance_values.append(distance)

        estimated_infection_rate = (1-p_values[numpy.argmin(distance_values)])
        
        return estimated_infection_rate, distance_values
    
    elif method == Infection_Rate_Method.RELATIVE_ECDF_AREA:
        
        p_values = numpy.linspace(0, 1, resolution+1)

        unique_counts = set(background_transcript_counts).union(set(signal_transcript_counts)).difference([0])
        unique_counts = sorted(list(unique_counts))

        distance_values = []

        distances_df_data = []

        for p_noise_estimate in p_values:

            distance = 0
            previous_count = 0
            
            previous_p_noise = background_transcript_counts[background_transcript_counts == 0].shape[0] * \
                p_noise_estimate/background_transcript_counts.shape[0]
            previous_p_signal = signal_transcript_counts[signal_transcript_counts == 0].shape[0]/\
                signal_transcript_counts.shape[0]

            previous_p_signal = min(previous_p_signal, p_noise_estimate)

            x_values = []
            noise_y_values = []
            signal_y_values = []

            last_above_index = 0

            for count_index, count in enumerate(unique_counts):

                x_values.append(count)

                # For this transcript count, calculate the expected probability of seeing it, assuming our debris percentage
                p_noise = background_transcript_counts[background_transcript_counts <= count].shape[0] * \
                    p_noise_estimate/background_transcript_counts.shape[0]
                noise_y_values.append(p_noise)

                # Similarly, calculate the probability of this transcript count in the debris + cell data
                p_signal = signal_transcript_counts[signal_transcript_counts <= count].shape[0]/\
                    signal_transcript_counts.shape[0]
                signal_y_values.append(p_signal)

                p_signal = min(p_signal, p_noise_estimate)

                points = []
                area_above = 0
                area_below = 0

                if previous_p_noise > previous_p_signal:
                    if p_noise > p_signal:
                        points.append([previous_count, previous_p_signal])
                        points.append([previous_count, previous_p_noise])
                        points.append([count, p_noise])
                        points.append([count, p_signal])
                        area_above = shoelace_formula(points)
                    else: # p_signal > p_noise
                        points.append([previous_count, previous_p_signal])
                        points.append([previous_count, previous_p_noise])
                        points.append([count, p_signal])
                        points.append([count, p_noise])
                        left_length = previous_p_noise - previous_p_signal
                        right_length = p_signal - p_noise
                        if left_length != 0 or right_length != 0:
                            height = count - previous_count
                            multiplier = height / (2*(left_length + right_length))
                            area_above = (left_length**2)*multiplier
                            area_below = (right_length**2)*multiplier
                else: # previous_p_signal > previous_p_noise
                    if p_signal > p_noise:
                        points.append([previous_count, previous_p_noise])
                        points.append([previous_count, previous_p_signal])
                        points.append([count, p_signal])
                        points.append([count, p_noise])
                        area_below = shoelace_formula(points)
                    else: # p_noise > p_signal
                        points.append([previous_count, previous_p_noise])
                        points.append([previous_count, previous_p_signal])
                        points.append([count, p_noise])
                        points.append([count, p_signal])
                        left_length = previous_p_signal - previous_p_noise
                        right_length = p_noise - p_signal

                        if left_length != 0 or right_length != 0:
                            height = count - previous_count
                            multiplier = height / (2*(left_length + right_length))
                            area_below = (left_length**2)*multiplier
                            area_above = (right_length**2)*multiplier

                distance += area_below + area_above

                previous_count = count
                previous_p_noise = p_noise
                previous_p_signal = p_signal

            if p_noise_estimate != 0:
                distance = distance/(p_noise_estimate*max(unique_counts))
            else:
                distance = 1

            distance_values.append(distance)

        estimated_infection_rate = (1-p_values[numpy.argmin(distance_values)])
        
        return estimated_infection_rate, distance_values
    
    elif method == Infection_Rate_Method.PAIRWISE_COMPARISONS:
        
        num_above = 0
        num_below = 0
        
        for count in signal_transcript_counts:
            num_above += (background_transcript_counts < count).sum()
            num_below += (background_transcript_counts > count).sum()
        
        estimated_infection_rate = (num_above-num_below)/(num_above+num_below)
        
        return estimated_infection_rate
    
    elif method == Infection_Rate_Method.NBINOM_P_NONZERO:
        
        signal_dist = get_neg_binomial_diffxpy(signal_transcript_counts)

        p_not_zero = 1 - signal_dist.pmf(0)

        return p_not_zero
    
    elif method == Infection_Rate_Method.NBINOM_FIXED_P:
        
        background_dist = get_neg_binomial_diffxpy(background_transcript_counts)
        
        signal_transcript_counts = numpy.round(signal_transcript_counts)

        signal_plus_background_params = scipy.optimize.minimize(
            neg_log_likelihood_nbinom_fixed_p,
            [1],
            args=(signal_transcript_counts, background_dist.args[1]),
            method="Nelder-Mead",
            options={
                "maxiter": 5000
            }
        )
        
        background_r = background_dist.args[0]
        signal_plus_background_r = numpy.exp(signal_plus_background_params["x"][0])
        p = background_dist.args[1]

        estimated_infection_rate = 1 - scipy.stats.nbinom(signal_plus_background_r - background_r, p).pmf(0)
        
        return estimated_infection_rate