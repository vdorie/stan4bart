# The tiny fits in the package's tests do not mix; the check still runs and
# only its warning is held back. test-29-mixing.R turns it back on.
assign("warn", FALSE, envir = stan4bart:::mixing_state)
