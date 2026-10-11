# B6 ch2-rehunt: D_riccati API fix + D2_DEEP_IM supersession note

Datum: 2026-10-11 (Abend). Autor: Bingsi (Hermes Agent). Baum: fix nach
`/home/error/physics/clones/SSZ-Transport-Bridge` gepatcht und in den
Hauptbaum transplantiert; Backups unter
`/root/deleted-backups/20261011-ch2-rehunt/`.

## 1. Befund: `D_riccati` war seit dem Refactor unaufrufbar

Refactor `46253bc -> 517ef2a` veraenderte die Outgoing-Selektion
("deterministic outgoing selection by |k - omega|" -> "3 largest Re lam")
und machte `Z_out(omega, r_eval, ...)` zwingend. `D_riccati(**kw)`
reichte `kw` aber an BEIDE Kandidaten weiter: `Z_in` nimmt kein
`r_eval` -> TypeError bei jedem Aufruf, vom try/except in
`sigma_min_riccati` verschluckt (Penalty-Rueckgabe 1e6).

Konsequenz: der 8.7-h Nelder-Mead-Polish-Lauf vom 11.10. (proc_c2d9874d65b8,
5 Starts am D2_DEEP_IM-argmin) war ein Totalausfall — alle 5 Starts
kehrten sofort mit Penalty zurueck, ohne eine einzige echte Funktionseval.
Ebenso veraltet: `tools/run_b6_step1_blind_scan.py` (Aufruf ohne r_eval;
Lauf Rueckgabe nur durch den Legacy-Default-Zweig vor dem Fix).

## 2. Fix

`D_riccati(omega, r_eval=None, **kw)`: explizite Trennung — `Z_in`
ohne `r_eval`, `Z_out` mit `r_eval` (Default `r_match = 30`, Konvention
der zertifizierten Controls). Verifikation: Endwerte stabil unter
rtol 1e-7 == 1e-9; 54/54 Tests gruen nach Transplant in den Hauptbaum.

## 3. Supersession-Kette: D2_DEEP_IM-Kandidat ist tot

Fokussierter Re-Scan (42 Punkte um den alten argmin (0.05, -1.4273),
aktuelle Maschinerie, 2026-10-11, Log
`/root/.hermes/cache/scratch/ch2_rehunt_local_scan.log`):

- sigma_min bleibt ueberall O(1): min 0.823 bei (0.10, -0.90);
  1.24 am alten argmin.
- r_eval-Sensitivitaet 30/45/60: 0.823 -> 0.821 -> 0.819 (regulaer,
  keine Divergenz Richtung 0).
- Der alte Artefaktwert sigma = 2.56e-3 (Blind-Scan 02:29 Uhr) stammt
  aus der VOR dem Refactor verwendeten Selektion und ist ueber den
  Refactor hinweg nicht vergleichbar.

Wichtig zur Einordnung: auch der zertifizierte ch2-Pilot-Pol
(0.2264259 - 0.0845509j) liefert unter sigma_min_riccati nur
sigma ≈ 0.59. Das WIDERSPRICHT der Zertifizierung nicht: sein eigenes
Artefakt (`SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json`) dokumentiert
`v1_riccati_3dof: SUPERSEDED (dichotomy collapse, wrong BVP measured)`.
Die Zertifizierung lief ueber den zweiseitigen vorwaerts-stabilen
QR-Evans-Wronskian-Piloten (Floors ~1e-21, r_match-Unabhaengigkeit
<1e-17, Mesh/rtol-Stabilitaet, Negative Controls).

## 4. Konsequenz

1. Der D2_DEEP_IM deep-IM-Kandidat (2.56e-3) ist unter der aktuellen
   Selektion NICHT reproduzierbar — und gehoerte ohnehin zur
   superseded Riccati-Linie. `SSZ_B6_STEP1_BLIND_SCAN_V1.json` bleibt
   als historisches Artefakt (vor-Refactor-Selektion) stehen; diese
   Datei dokumentiert die Supersession.
2. B6 Step 1 (physical pole) bleibt OPEN. Weitere Poljagd muss ueber
   die QR-Evans-Wronskian-Maschinerie laufen, nicht ueber
   sigma_min_riccati.
3. Der API-Fix ist Voraussetzung fuer jeden kuenftigen Riccati-gestuetzten
   Lauf; ohne ihn ist jede Penalty-Rueckgabe ein Maskierungsartefakt.
