function Q = general_WS(TM1, TM2, dx)

Q = -1i * TM2 * (TM2 - TM1)' / dx;

end