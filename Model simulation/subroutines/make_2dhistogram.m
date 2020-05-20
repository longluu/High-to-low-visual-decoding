function [est_hist_cond1, est_hist_cond2, est_hist_cond3_1, est_hist_cond3_2, est_hist_cond3_3] = ...
        make_2dhistogram(estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, smooth_sampling, edge_bin, h_filter)

% Make 2d histogram
if ~isempty(estimate_cond1)
    est_hist_cond1 = hist3(flipud(estimate_cond1)', 'Edges', {edge_bin, edge_bin});
    est_hist_cond2 = hist3(flipud(estimate_cond2)', 'Edges', {edge_bin, edge_bin});
else
    est_hist_cond1 = [];
    est_hist_cond2 = [];
end
if ~isempty(estimate_cond3_1)
    est_hist_cond3_1 = hist3(flipud(estimate_cond3_1)', 'Edges', {edge_bin, edge_bin});
    est_hist_cond3_2 = hist3(flipud(estimate_cond3_2)', 'Edges', {edge_bin, edge_bin});
    est_hist_cond3_3 = hist3(flipud(estimate_cond3_3)', 'Edges', {edge_bin, edge_bin});
end

% Smooth the histogram
if smooth_sampling
    if ~isempty(estimate_cond1)
        est_hist_cond1 = imfilter(est_hist_cond1, h_filter);
        est_hist_cond2 = imfilter(est_hist_cond2, h_filter);
    end
    if ~isempty(estimate_cond3_1)
        est_hist_cond3_1 = imfilter(est_hist_cond3_1, h_filter);
        est_hist_cond3_2 = imfilter(est_hist_cond3_2, h_filter);
        est_hist_cond3_3 = imfilter(est_hist_cond3_3, h_filter);
    end
end


