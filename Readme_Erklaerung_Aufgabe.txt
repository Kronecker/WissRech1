Anmerkung:
Das Programm wurde mit Clion und CMake erstellt. Das Projekt enthält auch die vorherigen Programmieraufgaben und kann unter https://github.com/Kronecker/WissRech1 komplett heruntergeladen werden.
git clone https://github.com/Kronecker/WissRech1.git

5 a)
Die Matrix ist dann diagonal dominant, wenn gilt:
v0 <= sqrt(5)/h  
(<= bedeutet kleiner gleich)
Detaillierte Erläuterung dazu befinden sich in "Erklaerung_Aufgabe_5a.jpg"

5 b) 
Das Programm erzeugt folgenden Output auf der Konsole: 

n       dd      v0      itera   time    Algo    Diff
25      y       50      67      0.5008  Jacobi  Central
25      y       50      21      0.501   Gauss   Central
25      n       100     197     1.5033  Jacobi  Central
25      n       100     85      1.0033  Gauss   Central
25      n       200     753     4.0099  Jacobi  Central
25      n       200     9999    53.6435 Gauss   Central
50      y       50      307     6.6222  Jacobi  Central
50      y       50      131     3.0383  Gauss   Central
50      y       100     131     3.0426  Jacobi  Central
50      y       100     35      1.0021  Gauss   Central
50      n       200     9999    220.587 Jacobi  Central
50      n       200     9999    220.085 Gauss   Central
100     y       50      1273    101.688 Jacobi  Central
100     y       50      603     56.5684 Gauss   Central
100     y       100     569     45.9377 Jacobi  Central
100     y       100     227     21.5569 Gauss   Central
100     y       200     253     21.0181 Jacobi  Central
100     y       200     59      6.3906  Gauss   Central
25      y       50      133     1.0023  Jacobi  Upwind
25      y       50      59      1.0304  Gauss   Upwind
25      y       100     85      1.0033  Jacobi  Upwind
25      y       100     31      0.5017  Gauss   Upwind
25      y       200     61      0.5013  Jacobi  Upwind
25      y       200     19      0.5018  Gauss   Upwind
50      y       50      433     9.4756  Jacobi  Upwind
50      y       50      201     5.4159  Gauss   Upwind
50      y       100     249     5.0129  Jacobi  Upwind
50      y       100     99      2.5069  Gauss   Upwind
50      y       200     161     3.5094  Jacobi  Upwind
50      y       200     51      1.0033  Gauss   Upwind
100     n       50      1525    120.861 Jacobi  Upwind
100     n       50      741     68.1819 Gauss   Upwind
100     y       100     797     63.6691 Jacobi  Upwind
100     y       100     349     36.0957 Gauss   Upwind
100     y       200     463     37.0988 Jacobi  Upwind
100     y       200     171     16.7236 Gauss   Upwind

Hierbei steht dd für diagonal dominant, itera sind die Anzahl der Iterationen (limitiert auf 9999), time ist die gebrauchte Zeit in s, Algo gibt an ob die Jacobi Iteration oder das Gauss Seidel Verfahren zum Lösen des LGS genutzt wurde, Diff gibt an ob zentrale oder Upwind Differenzen.
Des Weiteren werden alle Ergebnisse in Dateien der Form "T5GaussCentral_n=0025_vo=050.dat" wobei der Name selbsterklärend ist. Die Dateien können mit gnuplots "set pm3d" und "splot" eingelesen werden. Ein vorgefertigtes Skript für eine Windows basierte Gnuplot Installation wird ebenfalls erzeugt.
<to be continued>

5 c)
Während die Anwendung der zentralen Differenzen bei n=100 und v0=200 eine nicht diagonal dominante Matrix erzeugt, erzeugen die Upwind Differenzen durchaus eine diagonal dominante Matrix. Ersteres lässt sich mit Jacobi oder Gauss Seidel entsprechend nicht sinnvoll lösen, während letzteres eine brauchbare Lösung erzeugt.



