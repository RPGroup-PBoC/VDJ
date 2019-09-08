// Define useful variables
var mut = drop.value;
filter.group = mut

// Set filters on view data
seq_view.filters[0] = filter 

// Update displayed view
seqs.data.view = seq_view

// Push change
seqs.change.emit()
