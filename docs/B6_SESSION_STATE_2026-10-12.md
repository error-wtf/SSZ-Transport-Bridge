# B6 Session-State (Snapshot 12.10.2026, für Session-Fortsetzung nach /new)

## Laufender Prozess
- PID 1803215 (100% CPU, gestartet 11.10.), `tools/run_b6_coupled_evans_scan.py`
  in /home/error/physics/SSZ-Transport-Bridge (.venv), v2-Scan.
- Log: /root/.hermes/cache/scratch/b6_coupled_scan_v2.log (SCAN_EXIT wird angehängt).
- WICHTIG: Prozess nutzt den beim Start geladenen Code; Grid-Persistenz-Änderung
  wirkt erst beim NÄCHSTEN Start. Prozess NICHT abwürgen, natürlich auslaufen lassen.

## Ergebnisse Stand 12.10. ~03:00
- Grid 294/294 Punkte erfolgreich. Pilot-Anker D=3.027e-09.
- Grid-Minimum: 8.075e-17.
- Refined Kandidaten (aus Log):
  1. w = 0.0518750000 + 0.0285000000j, D = 5.501e-18  (kleiner als Grid-Min)
  2. w = 0.0499768066 + 0.2101666260j, D = 7.871e-17
- Beide im Bereich numerischer Präzisionsgrenzen — Kleinheit allein beweist KEINE Nullstelle.

## Offene Prüfpunkte (nächste Session)
1. ZEITKONVENTION im Solver klären: e^{-iωt} vs e^{+iωt}. Bei e^{-iωt} wäre
   Im ω > 0 ein WACHSTUMS-/Instabilitätsmodus, kein gedämpfter QNM-Pol.
   Grep-Hinweis: tools/run_b6_coupled_evans_scan.py Zeile ~6: "half w-plane
   (scan convention: poles reported as conjugates, exact ...)" — genau lesen!
2. Argumentprinzip / Nullstellenzählung auf kleiner Kontur um Kandidat 1
   (komplexes D(ω), nicht nur |D|).
3. Gates: r_match 25/30/35, r_seed 50/55/59.5, Konjugationsbeziehung,
   Negativkontrollen (Codezeilen ~188-198) — Werte aus finale
   COMPLETE-Artefakt auswerten.
4. Finale Artefakte landen in /home/error/physics/SSZ-Transport-Bridge/artifacts/
   (Verzeichnis existiert, mtime 12.10. 02:54 — schon reinsehen).
5. tests/ Phase & Observable Contract V1: 74 Tests PASS, Ruff sauber (laut
   vorheriger Session). Physische Quelle: Quellecouplung aus gemeinsamer
   Wirkung weiterhin OFFEN.

## Bewertungshalt (mit ChatGPT abgestimmt)
- Entscheidend ist NICHT 1e-18, sondern: robuste, verifizierte komplexe
  Evans-Nullstelle + korrekte physikalische Klassifikation (Gate).
- Validierter Nullpunkt > bloße Residuen-Suche.

## Addendum 12.10., 06:20 Uhr (Bingsi, nach Carmens Klassifikations-Fragen)

**Zeitkonvention verifiziert** (independent_qnm_3dof_v1.py Docstring + stable_evans_v3.py:201):
- ODE in t-time convention, omega = omega_t.
- Scan-Konvention: Nullstellen werden bei Im w > 0 gesucht und als
  KONJUGIERTE Pole gemeldet (`pole_qnm_sign` mit Im<0 im Artefakt).
- => Kandidat w = 0.051875 + 0.0285i (Scan-Form) entspricht einem
  QNM-Pol 0.051875 - 0.0285i => e^{-iωt}: DAMPFT. Kein
  Instabilitätsmodus per Konvention. (Sachliche Prüfung der
  Branch-Wahl bleibt Teil der Gate-Auswertung beim COMPLETE.)

**Argumentprinzip V1 (02:54) ist UNGÜLTIG — Probe hatte Bug:**
- Closure-Term `angle(v_last*conj(v0))` halbiert positive Windungen
  ( analytischer Check: W=w-w0 gab N=+0.5 statt +1 ) und kann
  Vorzeichen kippen (Kontrolle: N=-1.33 statt 0).
- Positivkontrolle Pilot-Pol scheiterte (N=0.139 < 1) => v1-Verdict
  NOT_CONFIRMED ist wertlos, MIPS der Methode, nicht des Kandidaten.
- Fix: `angle(v0*conj(v_last))`; Sanitätstest an analytischen
  Funktionen 5/5 PASS (+1/0/-1/+2/0).
- v2-Lauf mit per_side=32 (statt 8): PID 2457354,
  Log /root/.hermes/cache/scratch/b6_argprinciple_v2.log, ETA ~2.4h.
  Konturen: candidate h=0.004/0.008/0.015, pilot h=0.008, offset h=0.008.

**Haupt-Scan-Status** (PID 1803215, gestartet 11.10. ~20:45):
- LEBT noch (99.9% CPU, Nelder-Mead in Kandidat 3/5).
- 3 verfeinerte Kandidaten bisher:
  1. w=0.051875+0.0285i        D=5.501e-18
  2. w=0.0499768066+0.2101666i D=7.871e-17
  3. w=0.0521972656+0.1479492i D=2.214e-17
- WICHTIG: Der laufende Prozess nutzt noch den Code VOR der Grid-
  Persistenz-Änderung (Linos Regel 12.10.) — daher kein
  SSZ_B6_STEP1_COUPLED_EVANS_GRID_V1.npz und keine "grid persisted"-
  Zeile. Bei Crash wären die 294 Grid-Punkte verloren (nur die
  refined-Kandidaten + top20 überleben im Log/Artefakt). Persistenz
  wird beim NÄCHSTEN Lauf wirksam.
- Danach: Gates pro Kandidat (r_match 25/30/35, r_seed 50/55/59.5,
  Konjugation, Negativkontrolle) => COMPLETE-Artefakt
  artifacts/SSZ_B6_STEP1_COUPLED_EVANS_SCAN_V1.json.

**Backups:** /root/deleted-backups/20261012-061150-b6-argprinciple-fix/
(Probe v1 + v1-Artefakt).

## Addendum 2, 12.10. 07:05 — Delta-Skalierung + SVD-Kontrolle (Carmens Punkte 2+3)

Artefakt: artifacts/SSZ_B6_STEP1_DELTA_SCALING_V1.json (45,614 B,
SHA-NPZ im File). 4 Zentren x (20 Ladder-Punkte in 2 Richtungen +
D_center + 4 rtol-Sweeps mit vollem sval-Spektrum), identische Physik
wie Scan (r_match=30, r_seed=55, 40 Segmente, rtol_ref=1e-9).

### Kernbefunde

1. KEINER der 3 Kandidaten zeigt Nullstellen-Struktur:
   - dir0/dir1 Steigungen m: cand1 0.10/0.01, cand2 0.02/-0.08,
     cand3 0.16/0.04 — erwartet waere m~1 fuer eine einfache NS.
   - |D| bleibt auf der GESAMTEN Leiter (delta 1e-2 .. 3e-8) auf
     1e-16..1e-14 skaliert (SVD-Roundoff-Floor ~ eps*s1 ~ 1e-16).
     Kein Anstieg nach aussen => kein Dip, nur Plateau.
2. rtol-Sweep 1e-9 -> 1e-10 -> 1e-11: D bit-identisch an allen
   Zentren (4 signifikante Stellen). ODE-Diskretisierung ist NICHT
   der dominante Fehler; der Floor ist SVD-/Subspace-Roundoff,
   rtol-unabhaengig. Carmens Punkt 4 (echte Hohe Praezision) ist
   damit die saubere Entscheidung, NICHT noch mehr rtol-Gedreh.
3. POSITIVKONTROLLE Pilot-Pol: D_center=3.027e-09, Umgebung
   1e-8..1e-7, ebenfalls kein m~1 ueber die Leiter. Die gekoppelte
   3-Kanal-Maschine loest den zertifizierten Pol nur bis ~1e-9 auf
   (vs. 1e-21 im 1x1-Pilot). => Die Maschine hat eine effektive
   Aufloesungsgrenze ~1e-9; ALLES darunter (alle 3 Kandidaten!) liegt
   unterhalb der demonstrierten Verlaesslichkeit der Maschine.
4. sval-Spektrum: [1.41, 1.41, 1.00, 1.00, ~0.03, Floor] an allen
   Kandidaten; Pilot: 5. sval 0.212. Die Luecke s5->s6 ist 13+ OD.

### Vorlaeufige Klassifikation (wartet auf ArgPrinciple v2 + Gates)

Die 3 Kandidaten sehen aus wie numerische Plateaus (Noise-Valleys),
nicht wie Resonanznullstellen: Grid-Min 8.075e-17 und refined
5.5e-18/2.2e-17/7.9e-17 liegen ALLE am SVD-Roundoff-Floor, weit
unter der maschinellen Aufloesung ~1e-9 (Pilot-Befund).

Konsistenz-Check mit Preregisterung: Die Negativkontrolle
(|D| bei w+0.01 muss DEUTLICH ueber Floor liegen) scheitert fuer
alle 3 (Ladder delta=1e-2: 6.0e-15 / 2.7e-13 / 5.4e-14 — gleiche
Groessenordnung wie Center). => Die vorregistrierten Gates sollten
die Kandidaten als REJECTED klassifizieren, wenn der Haupt-Scan
dort ankommt. Das ist das System working as designed.

Offen:
- ArgPrinciple v2 (PID 2457354) laeuft, ETA ~08:30. Erwartung bei
  Plateau: N~0 (keine Windung). Faellt N=1 an einem Kandidaten aus,
  waere das ein echter Widerspruch -> untersuchen.
- Haupt-Scan (PID 1803215): NM Kandidat 3/5, danach Gates.
- Evtl. Follow-up (Carmen Punkt 4): Propagation+Matching wirklich in
  mpmath (30-50 digits) neu implementieren und die Pilot-Region +
  eine Kandidat-Region nachrechnen. Das wuerde die wahre Landschaft
  zeigen (double-Roundoff raus). Aufwand: neu, noch nicht geplant.

## Addendum 3, 12.10. 07:35 — SCAN COMPLETE + ArgPrinciple v2 + finale Methoden-Bilanz

### 1. Haupt-Scan: COMPLETE (05:55 UTC, 11.1 h, 294 grid / 317 evals / 0 errors)
FIVE Kandidaten (nicht 3), alle classification GATED_POLE_CANDIDATE,
alle gate_pass=true: w_scan = 0.051875+0.0285i (5.5e-18),
0.049977+0.21017i (7.9e-17), 0.052197+0.14795i (2.2e-17),
0.050742+0.32742i (5.0e-17), 0.050624+0.38435i (5.9e-17).
Realteile ALLE ~0.05; Imaginaerteile 0.03..0.38 (Modenleiter-Verdacht).

### 2. ArgPrinciple v2 (per_side=32, gefixte Closure): NOT_CONFIRMED —
und zwar BEWEISBAR aus Methodengruenden, nicht wegen der Kandidaten:
- Pilot-Positivkontrolle (zertifizierter Pol, |W|~1e-10..1e-7, also
  7+ OD ueber Floor): N=0.0 — echte Null wuerde N=1 winden.
- Offset-Kontrolle: N=-3.0 (muss 0 sein).
- Kandidat: N = -0.0 / +3.0 / +1.0 bei h=0.004/0.008/0.015 —
  nicht mal konturkonsistent => wertlos.
ROOT CAUSE (Diagnose tools/b6_pilot_phase_diag.py): W=det([Qi|Qe])
nach Spalten-Normalisierung + QR-Sign-Fix ist KEIN analytisches W in
omega: die Phase springt zwischen benachbarten Konturpunkten um
~1-4 rad (arg-Sequence im Log), |W| bleibt glatt. Argumentprinzip
braucht eine analytische (oder stetig analytisch fortgesetzte)
Evans-Funktion; die normalisierte det-Version ist stueckweise definiert
und Phasenspruenge sind Normalisierungsartefakte. => Windungsprobe
braucht Neuimplementierung (rohe det ohne Normalisierung ODER
Compound-Matrix/Exterior-Product-Evans-Funktion mit Phasen-Tracking).
Probe in aktueller Form RETIRED.

### 3. Gate-Design-Befund (Scan): Die vorregistrierten Gates sind
skalen-relativ und werden am Roundoff-Floor vakuum-passiert:
d_ref < FLOOR_TARGET(1e-12) - von Noise voellig erfuellt; negctrl > 100*max(d_ref,1e-18):
bei d_ref=5.5e-18 reicht 5.5e-16, Negativkontrolle ist 1.3e-14
(=2400x) -> PASS, obwohl 1.3e-14 selbst nur Floor+1 OD ist.
Gegenprobe: der zertifizierte Pilot wuerde den negctrl-Gate FAILEN
(D(0.01 weg) nur ~30x D_center). => GATED_POLE_CANDIDATE ist hier
keine physikalische Klassifikation, sondern ein Threshold-Artefakt.

### 4. Finale Methoden-Triangulation (12.10., Stand 07:35)
- delta-Skalierung: KEIN Dip an allen 5 Kandidaten-Umgebungen; alle
  |D| = 1e-16..1e-18 = SVD-Roundoff-Floor; Machine-Loesungsgrenze
  (am Pilot demonstriert) ~3e-9. Kandidaten 6-8 OD DARUNTER.
- rtol 1e-9/1e-10/1e-11: bit-identisch => Floor ist SVD-Roundoff,
  nicht Diskretisierung.
- ArgPrinciple: in dieser Implementierung unbrauchbar (Phasen-Normal-
  isierungsartefakte), bewiesen an der Positivkontrolle.
- Gates: vakuum-passiert (s.o.).
SCHLUSSFOLGERUNG (vorlaeufig, bis mpmath-Nachrechnung): Die 5
Kandidaten sind mit hoher Wahrscheinlichkeit numerische Noise-Valleys
am SVD-Floor, NICHT validierte Evans-Nullstellen. Der einzige
Maschinen-bestaetigte Pol bleibt der Pilot (D~3e-9).

### 5. Offen / Naechste Schritte (Entscheidung bei Lino/Carmen)
a) mpmath-Propagator (30-50 digits) fuer Pilot-Region + eine
   Kandidat-Region: zeigt, ob unter dem Double-Floor ueberhaupt
   analytische Landschaft liegt (Carmens Punkt 4, jetzt begruendet).
b) ArgPrinciple-Neubau mit roher det / Compound-Matrix.
c) Gate-Rekalibrierung: absolute Machine-Floor-Schwelle (z.B.
   d_ref < 1e-12 UND negctrl > 1000x Floor UND Dip-Tiefe >= 2 OD
   gegen Umgebung) statt skalen-relativer Schwellen.

---

## Addendum 4 (12.10. ~13:20) — mpmath-Kontrollsolver: Drei-Fälle-Test KOMPLETT

**Solver**: `tools/b6_mpmath_evans_v2.py` (rohe det von M_match, holomorph
by construction — kein QR in der Funktion; Spalten-normalisiertes σ_min als
skaleninvarianter Indikator). Verifikation: ∏σ = |det| auf 15+ Stellen,
RK4-Schrittverdopplung driftet σ_min nur ~15%, r_match-Sensitivität klein,
incoming-µ-Auswahl im Pilot-Gebiet kontinuierlich (keine Basiswechsel).

**Ergebnis des Drei-Fälle-Tests (dps 30/50/80):**

| Fall | σ_min | Dip bei δ→0 | Urteil |
|---|---|---|---|
| Pilot 0.22643+0.08455i | 8.6e-9 (roh), 7.2e-9 (reg) | NEIN | **KEINE Nullstelle** |
| Kandidat 1 (Re≈0.05) | 1.4e-5 | NEIN (nicht-monoton) | Keine Nullstelle |
| Negativkontrolle | 6.5e-6 | NEIN | Sauber — Kandidat nicht unterscheidbar |

5×5-Regions-Grid (Re 0.15–0.35, Im 0.03–0.23): σ_min_reg fällt monoton zur
Ecke (0.35, 0.23) → Minimum 1.5e-11, **kein lokaler Dip am Pilot**.

**Zentrale Kehrtwende**: Der "zertifizierte Pilot-Pol" ist **keine Nullstelle
von E**. Die Double-Präzision hörte bei 3e-9 auf (Floor-Illusion); mpmath
bewegt sich Glatt DURCH 3e-9 hindurch ohne jeden Pol. Hohe Präzision hat die
Frage nicht "schärfer beantwortet", sondern die frühere Bestätigung
widerlegt. Die fünf Re≈0.05-Kandidaten bleiben Rauschen (Rosette
dps-unabhängig, Niveau = Negativkontrolle).

**Validitätsgrenze**: jenseits t≈0.3 auf der Diagonale kollabiert der
Frobenius-Seed auf Rang-1-Dominanz (log_scale_int 7→186; σ₃ fällt mit,
2 Null-svals bei t=1.0) — der weitere Abfall dort ist Konstruktions-Artefakt.

**Offen**: Wo (ob) E im Gültigkeitsbereich wirklich verschwindet;
Gate-Rekalibrierung; ggf. Argumentprinzip auf der rohen holomorphen det mit
Phasen-Tracking als Nullstellen-Freiheits-Zertifikat.

**Artefakte**: `SSZ_B6_STEP1_MPMATH_EVANS_{PILOT,CANDIDATE,NEGATIVE}_DPS{30,50,80}_R1.json`
+ `SSZ_B6_STEP1_MPMATH_CONTROLSOLVER_VERDICT_V1.json` (Gesamturteil).

---

## Addendum 5 (12.10.26, ~15:00) — Carmen-Priorisierung umgesetzt: Supersession + Gate-Rekalibrierung + Solver-Validierung

**1. Pilot-Zertifikat formell SUPERSEDET** (`SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json` in-place markiert + `SSZ_B6_PILOT_SUPERSESSION_V1`-Block im gleichen File):
- Status CANDIDATE_POLE_CERTIFIED_NUMERICAL → supersedet durch mpmath-Kontrollsolver (`7215aec`).
- Pilot 0.22643+0.08455i ist KEINE gültige Referenz/Validierungs-Anker mehr für B6. Die Reproduktion im Coupled-Scan zählt nur noch als Machinery-Cross-Check.
- Scope: nur der Pole-Claim und die Gates; FK/Sagnac/Transport/lokale Spektralresultate ausdrücklich NICHT berührt.

**2. Gates rekalkibriert** (`SSZ_B6_GATE_RECALIBRATION_V1` + Code-Änderung in `tools/run_b6_coupled_evans_scan.py`):
- Alt (vakuum, alle 5 Kandidaten vom 12.10. Rauschen): absolut |D|<1e-12 + relative Negativkontrolle `100×max(d_ref,1e-18)`.
- Neu: gemessener lokaler Floor (Rosette, 4 Offsets, zweitkleinster), Dip ≥ 2 OD, sichtbare Skalierung 50×/Dekade, absolute Negativkontrolle gegen den GEMESSENEN Floor, `gate_version=RECALIBRATED_V1`. Finale Bestätigung NUR via mpmath-Kontrollsolver (Rosette + dps-Leiter).

**3. Solver gegen analytisch bekannte Nullstellen validiert** (`tools/b6_solver_validation_analytic_v1.py` → `SSZ_B6_STEP1_SOLVER_VALIDATION_V1.json`): **VALIDATION_PASSED** (exit 0)
- Halblinien-Kasten (V0=-2.5, a=1), gleiche Konstruktionskette wie Kontrollsolver (RK4-Propagation, Matching an r=a, rohe holomorphe det, svals via eig(M^H M)).
- Einfache Nullstelle (N1-Referenzwurzel, 100-stellig): Dip 9.9 OD, Slope 1.000, R²=1.0 → PASS.
- Doppelpol-Exceptional-Point: analytisch reduziert auf w=-i/a mit q·tan(q·a) = -1/a (a0=1.8559…, w0=-0.538798…i; Residuen F, dF/dw ~ 1e-51): Dip 6.8 OD, Slope 2.02, R²=0.99994 → PASS mit deklarierter Ordnung 2.
- Negativkontrolle: REJECTED (dip −0.0 OD). Retro-Test: alter Pilot durch neue Gates → REJECTED (dip −0.11 OD).
- Konstruktions-Invarianten: ∏svals=|det| Fehler 2.2e-24; Referenzwurzel im.closed-form-Residuum 1.0e-16.
- Fix-Historie im Lauf: (a) Regressions-Bug (Zentrum mit δ=0 im Fit) — Zentrum aus G2-Fit entfernt; (b) erste Doppelnull-Algebra falsch (residual_dF_dw=4.8 entlarvt) → korrekte Reduktion w=-i/a.

**4. Nächster Schritt (offen)**: Physik-SSZ-Suche mit rekalibrierten Gates + holomorpher Evans-Funktion; optional ArgPrinciple v3 auf roher det mit Phasen-Tracking (QR-Version war der Bug, nicht das Prinzip).
