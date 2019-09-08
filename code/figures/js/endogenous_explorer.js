// Define useful variables
var mut = sel.value;
filter.group = mut;


// Set filters on view data
seq_view.filters[0] = filter;
dwell_view.filters[0] = filter;
cut_view.filters[0] = filter;
post_view.filters[0] = filter;
// Update displayed view
seqs.data.view = seq_view;
dwell_hist.data.view = dwell_view;
cut_hist.data.view = cut_view;
post.data.view = post_view;

// Push change
seqs.change.emit();
dwell_hist.change.emit();
cut_hist.change.emit();
post.change.emit()
