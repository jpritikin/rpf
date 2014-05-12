set print thread-events off
dir src
set breakpoint pending on
b Rf_error
# Hints:
# Use Rf_PrintValue(SEXP) to pretty print R data from gdb
