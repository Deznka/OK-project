// main4.cpp - Zoptymalizowany algorytm mrówkowy ACO dla TSP
// Cel: jak najlepsze wyniki w max 3 minuty dla instancji do 1000 miast

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <random>

using namespace std;

// ---------------------------------------------------
// GLOBALNE ZMIENNE CZASOWE
// ---------------------------------------------------
chrono::steady_clock::time_point globalStartTime;
const int MAX_TIME_SECONDS = 175; // bezpieczny limit (< 3 min)

bool czasPrzekroczony() {
    auto teraz = chrono::steady_clock::now();
    auto sek = chrono::duration_cast<chrono::seconds>(teraz - globalStartTime).count();
    return sek >= MAX_TIME_SECONDS;
}

int pozostalyCzas() {
    auto teraz = chrono::steady_clock::now();
    auto sek = chrono::duration_cast<chrono::seconds>(teraz - globalStartTime).count();
    return max(0, MAX_TIME_SECONDS - (int)sek);
}

// ---------------------------------------------------
// KOSZT TRASY
// ---------------------------------------------------
double policzKosztTrasy(const vector<int>& trasa, const vector<vector<double>>& D) {
    double koszt = 0.0;
    int n = (int)trasa.size();
    for (int i = 0; i < n - 1; i++) {
        koszt += D[trasa[i]][trasa[i + 1]];
    }
    koszt += D[trasa[n-1]][trasa[0]];
    return koszt;
}

// ---------------------------------------------------
// SZYBKI 2-OPT Z LIMITEM ITERACJI
// ---------------------------------------------------
void twoOptFast(vector<int>& trasa, const vector<vector<double>>& D, int maxPasses = 100) {
    int n = (int)trasa.size();
    if (n < 4) return;

    for (int pass = 0; pass < maxPasses; pass++) {
        bool poprawa = false;

        for (int i = 1; i < n - 2 && !poprawa; i++) {
            for (int j = i + 1; j < n - 1; j++) {
                int a = trasa[i - 1];
                int b = trasa[i];
                int c = trasa[j];
                int d = trasa[(j + 1) % n];

                double przed = D[a][b] + D[c][d];
                double po = D[a][c] + D[b][d];

                if (po + 1e-9 < przed) {
                    reverse(trasa.begin() + i, trasa.begin() + j + 1);
                    poprawa = true;
                    break;
                }
            }
        }

        if (!poprawa) break;
    }
}

// ---------------------------------------------------
// PEŁNY 2-OPT (dla małych instancji)
// ---------------------------------------------------
void twoOptFull(vector<int>& trasa, const vector<vector<double>>& D) {
    int n = (int)trasa.size();
    if (n < 4) return;

    while (true) {
        bool poprawa = false;

        for (int i = 1; i < n - 2 && !poprawa; i++) {
            for (int j = i + 1; j < n - 1; j++) {
                int a = trasa[i - 1];
                int b = trasa[i];
                int c = trasa[j];
                int d = trasa[(j + 1) % n];

                double przed = D[a][b] + D[c][d];
                double po = D[a][c] + D[b][d];

                if (po + 1e-9 < przed) {
                    reverse(trasa.begin() + i, trasa.begin() + j + 1);
                    poprawa = true;
                    break;
                }
            }
        }

        if (!poprawa) break;
    }
}

// ---------------------------------------------------
// 3-OPT DLA MAŁYCH INSTANCJI (dodatkowa optymalizacja)
// ---------------------------------------------------
void threeOptSegment(vector<int>& trasa, const vector<vector<double>>& D, int maxIterations = 50) {
    int n = (int)trasa.size();
    if (n < 6) return;

    for (int iter = 0; iter < maxIterations; iter++) {
        bool improved = false;

        for (int i = 0; i < n - 4 && !improved; i++) {
            for (int j = i + 2; j < n - 2 && !improved; j++) {
                for (int k = j + 2; k < n; k++) {
                    // Sprawdź różne rekombinacje segmentów
                    int a = trasa[i], b = trasa[i+1];
                    int c = trasa[j], d = trasa[j+1];
                    int e = trasa[k], f = trasa[(k+1) % n];

                    double current = D[a][b] + D[c][d] + D[e][f];
                    double opt1 = D[a][c] + D[b][d] + D[e][f]; // 2-opt i-j
                    double opt2 = D[a][b] + D[c][e] + D[d][f]; // 2-opt j-k

                    if (opt1 + 1e-9 < current) {
                        reverse(trasa.begin() + i + 1, trasa.begin() + j + 1);
                        improved = true;
                        break;
                    }
                    if (opt2 + 1e-9 < current) {
                        reverse(trasa.begin() + j + 1, trasa.begin() + k + 1);
                        improved = true;
                        break;
                    }
                }
            }
        }
        if (!improved) break;
    }
}

// ---------------------------------------------------
// NEAREST NEIGHBOR - SZYBKA TRASA STARTOWA
// ---------------------------------------------------
vector<int> nearestNeighbor(int n, const vector<vector<double>>& D, int start = 1) {
    vector<bool> odw(n + 1, false);
    vector<int> trasa;
    trasa.reserve(n);

    trasa.push_back(start);
    odw[start] = true;

    int akt = start;
    for (int k = 1; k < n; k++) {
        int best = -1;
        double bestDist = numeric_limits<double>::infinity();
        for (int j = 1; j <= n; j++) {
            if (!odw[j] && D[akt][j] < bestDist) {
                bestDist = D[akt][j];
                best = j;
            }
        }
        if (best == -1) break;
        trasa.push_back(best);
        odw[best] = true;
        akt = best;
    }
    return trasa;
}

// ---------------------------------------------------
// STRUKTURY PARAMETRÓW
// ---------------------------------------------------
struct ACOParams {
    int liczbaMrowek;
    int iteracje;
    double alfa;
    double beta;
    double rho;
    double Q;
    int eliteFactor;
    int candidateSize;
    double tauMin;
    double tauMax;
    int opt2Passes;  // ile przejść 2-opt
};

struct ACOResult {
    vector<int> najlepszaTrasa;
    double najlepszyKoszt;
};

// ---------------------------------------------------
// PARAMETRY ADAPTACYJNE ZALEŻNE OD ROZMIARU PROBLEMU
// ---------------------------------------------------
ACOParams getAdaptiveParams(int n) {
    ACOParams p;

    // UNIWERSALNE PARAMETRY bazowane na testach
    // Kluczowe odkrycie: niska beta (1-2) + wysokie rho (0.5-0.8) + niskie Q (10-50)
    // daje lepsze wyniki niż klasyczne ustawienia (beta=5, rho=0.3, Q=100)

    if (n <= 60) {
        // Małe instancje (berlin52)
        p.liczbaMrowek = min(n, 50);
        p.iteracje = 600;
        p.alfa = 1.0;
        p.beta = 2.0;      // niska - więcej eksploracji
        p.rho = 0.6;       // wysokie parowanie - szybka adaptacja
        p.Q = 10.0;        // niskie Q
        p.eliteFactor = 3;
        p.candidateSize = min(20, n - 1);
        p.tauMin = 0.01;
        p.tauMax = 6.0;
        p.opt2Passes = 300;
    }
    else if (n <= 150) {
        // Średnie instancje (bier127)
        p.liczbaMrowek = min(n / 3, 40);
        p.iteracje = 500;
        p.alfa = 1.0;
        p.beta = 1.0;      // bardzo niska - sprawdzone na bier127
        p.rho = 0.7;       // wysokie parowanie
        p.Q = 10.0;        // niskie Q
        p.eliteFactor = 3;
        p.candidateSize = min(25, n - 1);
        p.tauMin = 0.01;
        p.tauMax = 5.0;
        p.opt2Passes = 200;
    }
    else if (n <= 300) {
        // Większe instancje (tsp250)
        p.liczbaMrowek = min(n / 4, 50);
        p.iteracje = 400;
        p.alfa = 1.0;
        p.beta = 1.5;
        p.rho = 0.65;
        p.Q = 15.0;
        p.eliteFactor = 4;
        p.candidateSize = min(30, n - 1);
        p.tauMin = 0.005;
        p.tauMax = 4.0;
        p.opt2Passes = 120;
    }
    else if (n <= 600) {
        // Duże instancje (tsp500)
        p.liczbaMrowek = min(n / 8, 60);
        p.iteracje = 300;
        p.alfa = 1.0;
        p.beta = 2.0;
        p.rho = 0.6;
        p.Q = 20.0;
        p.eliteFactor = 5;
        p.candidateSize = min(35, n - 1);
        p.tauMin = 0.001;
        p.tauMax = 3.0;
        p.opt2Passes = 60;
    }
    else {
        // Bardzo duże instancje (tsp1000)
        p.liczbaMrowek = min(n / 15, 60);
        p.iteracje = 200;
        p.alfa = 1.0;
        p.beta = 2.5;
        p.rho = 0.55;
        p.Q = 25.0;
        p.eliteFactor = 6;
        p.candidateSize = min(40, n - 1);
        p.tauMin = 0.0005;
        p.tauMax = 2.5;
        p.opt2Passes = 35;
    }

    // dla bier127
    // p.liczbaMrowek = min(n, 50);
    // p.iteracje = 500;
    // p.alfa = 1.0;
    // p.beta = 1;
    // p.rho = 0.7;
    // p.Q = 10.0;
    // p.eliteFactor = 5;
    // p.candidateSize = min(20, n - 1);
    // p.tauMin = 0.01;
    // p.tauMax = 10.0;
    // p.opt2Passes = 200;

    return p;
}

// ---------------------------------------------------
// GŁÓWNY ALGORYTM ACO
// ---------------------------------------------------
ACOResult runACO(const vector<vector<double>>& D, int n, const ACOParams& par,
                 unsigned int seed, bool verbose = false) {

    srand(seed);

    int liczbaMrowek = par.liczbaMrowek;
    int iteracje = par.iteracje;
    double alfa = par.alfa;
    double beta = par.beta;
    double rho = par.rho;
    double Q = par.Q;
    int eliteFactor = par.eliteFactor;
    int candidateSize = min(par.candidateSize, n - 1);
    double tauMin = par.tauMin;
    double tauMax = par.tauMax;
    int opt2Passes = par.opt2Passes;

    // Lista kandydatów (najbliższe miasta)
    vector<vector<int>> kandydaci(n + 1);
    for (int i = 1; i <= n; i++) {
        vector<pair<double,int>> tmp;
        tmp.reserve(n - 1);
        for (int j = 1; j <= n; j++) {
            if (i != j) tmp.push_back({ D[i][j], j });
        }
        sort(tmp.begin(), tmp.end());
        int kmax = min(candidateSize, (int)tmp.size());
        for (int k = 0; k < kmax; k++) {
            kandydaci[i].push_back(tmp[k].second);
        }
    }

    // Inicjalizacja feromonów
    double tau0 = 1.0;
    vector<vector<double>> feromony(n + 1, vector<double>(n + 1, tau0));

    // Atrakcyjność (1/odległość)
    vector<vector<double>> atrakcyjnosc(n + 1, vector<double>(n + 1, 0.0));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            if (i != j) atrakcyjnosc[i][j] = 1.0 / (D[i][j] + 0.0001);

    // Startowa trasa z Nearest Neighbor
    vector<int> najlepszaTrasa = nearestNeighbor(n, D, 1);
    twoOptFast(najlepszaTrasa, D, opt2Passes);
    double najlepszyKoszt = policzKosztTrasy(najlepszaTrasa, D);

    int bezPoprawy = 0;
    int maxBezPoprawy = max(20, iteracje / 3);

    // Pre-compute powers for speed
    vector<vector<double>> tauPow(n + 1, vector<double>(n + 1));
    vector<vector<double>> etaPow(n + 1, vector<double>(n + 1));

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            etaPow[i][j] = pow(atrakcyjnosc[i][j], beta);
        }
    }

    for (int iter = 0; iter < iteracje; iter++) {

        if (czasPrzekroczony()) {
            if (verbose) cout << "Stop: limit czasu globalny.\n";
            break;
        }

        // Update tau powers (po zmianie feromonów)
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                tauPow[i][j] = pow(feromony[i][j], alfa);
            }
        }

        vector<vector<int>> trasy(liczbaMrowek);
        vector<double> koszty(liczbaMrowek);
        bool poprawaWTejIteracji = false;

        // Każda mrówka buduje trasę
        for (int m = 0; m < liczbaMrowek; m++) {

            int start = (m % n) + 1;
            vector<int> trasa;
            trasa.reserve(n);
            vector<bool> odwiedzony(n + 1, false);

            trasa.push_back(start);
            odwiedzony[start] = true;
            int aktualne = start;

            for (int k = 1; k < n; k++) {

                double suma = 0.0;
                vector<pair<double, int>> opcje;
                opcje.reserve(candidateSize + 10);

                // Najpierw kandydaci
                for (int j : kandydaci[aktualne]) {
                    if (!odwiedzony[j]) {
                        double p = tauPow[aktualne][j] * etaPow[aktualne][j];
                        opcje.push_back({p, j});
                        suma += p;
                    }
                }

                // Jeśli wszyscy kandydaci odwiedzeni
                if (suma < 1e-10) {
                    for (int j = 1; j <= n; j++) {
                        if (!odwiedzony[j]) {
                            double p = tauPow[aktualne][j] * etaPow[aktualne][j];
                            opcje.push_back({p, j});
                            suma += p;
                        }
                    }
                }

                // Losowanie
                int wybrane = -1;
                if (suma > 1e-10) {
                    double los = ((double)rand() / RAND_MAX) * suma;
                    double akum = 0.0;
                    for (auto& op : opcje) {
                        akum += op.first;
                        if (akum >= los) {
                            wybrane = op.second;
                            break;
                        }
                    }
                }

                if (wybrane == -1 && !opcje.empty()) {
                    wybrane = opcje.back().second;
                }
                if (wybrane == -1) {
                    for (int j = 1; j <= n; j++)
                        if (!odwiedzony[j]) { wybrane = j; break; }
                }

                trasa.push_back(wybrane);
                odwiedzony[wybrane] = true;
                aktualne = wybrane;
            }

            // Lokalna optymalizacja 2-opt
            twoOptFast(trasa, D, opt2Passes);

            double koszt = policzKosztTrasy(trasa, D);
            trasy[m] = trasa;
            koszty[m] = koszt;

            if (koszt < najlepszyKoszt - 1e-9) {
                najlepszyKoszt = koszt;
                najlepszaTrasa = trasa;
                poprawaWTejIteracji = true;
            }
        }

        // Parowanie feromonów
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                feromony[i][j] *= (1.0 - rho);

        // Aktualizacja feromonów - wszystkie mrówki
        for (int m = 0; m < liczbaMrowek; m++) {
            double delta = Q / koszty[m];
            const auto& trasa = trasy[m];

            for (int i = 0; i < n - 1; i++) {
                int a = trasa[i];
                int b = trasa[i + 1];
                feromony[a][b] += delta;
                feromony[b][a] += delta;
            }
            feromony[trasa[n-1]][trasa[0]] += delta;
            feromony[trasa[0]][trasa[n-1]] += delta;
        }

        // Elitarne wzmocnienie najlepszej trasy
        double deltaElite = eliteFactor * (Q / najlepszyKoszt);
        for (int i = 0; i < n - 1; i++) {
            int a = najlepszaTrasa[i];
            int b = najlepszaTrasa[i + 1];
            feromony[a][b] += deltaElite;
            feromony[b][a] += deltaElite;
        }
        feromony[najlepszaTrasa[n-1]][najlepszaTrasa[0]] += deltaElite;
        feromony[najlepszaTrasa[0]][najlepszaTrasa[n-1]] += deltaElite;

        // Ograniczenie feromonów
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (feromony[i][j] < tauMin) feromony[i][j] = tauMin;
                if (feromony[i][j] > tauMax) feromony[i][j] = tauMax;
            }
        }

        if (verbose && (iter % 10 == 0 || poprawaWTejIteracji)) {
            cout << "Iter " << iter + 1 << " | Best: " << fixed << setprecision(2) << najlepszyKoszt << endl;
        }

        // Warunek stopu - brak poprawy
        if (poprawaWTejIteracji) bezPoprawy = 0;
        else bezPoprawy++;

        if (bezPoprawy >= maxBezPoprawy) {
            if (verbose) cout << "Stop: brak poprawy przez " << bezPoprawy << " iteracji.\n";
            break;
        }

        // Restart feromonów przy stagnacji
        if (bezPoprawy > 0 && bezPoprawy % (maxBezPoprawy / 2) == 0) {
            if (verbose) cout << "Restart feromonow (stagnacja).\n";
            for (int i = 1; i <= n; i++)
                for (int j = 1; j <= n; j++)
                    feromony[i][j] = tau0;
        }
    }

    // Finalna optymalizacja 3-opt dla małych instancji
    if (n <= 200) {
        threeOptSegment(najlepszaTrasa, D, 100);
        najlepszyKoszt = policzKosztTrasy(najlepszaTrasa, D);
    }

    ACOResult res;
    res.najlepszaTrasa = najlepszaTrasa;
    res.najlepszyKoszt = najlepszyKoszt;
    return res;
}

// ---------------------------------------------------
// ALGORYTM ZACHŁANNY (GREEDY) - wartość bazowa
// ---------------------------------------------------
double algorytmZachlanny(int n, const vector<vector<double>>& D, vector<int>& trasaOut) {
    // Znajdź najlepszą trasę startując z każdego miasta
    double najlepszyKoszt = numeric_limits<double>::infinity();
    vector<int> najlepszaTrasa;

    for (int startCity = 1; startCity <= n; startCity++) {
        vector<int> trasa = nearestNeighbor(n, D, startCity);
        double koszt = policzKosztTrasy(trasa, D);
        if (koszt < najlepszyKoszt) {
            najlepszyKoszt = koszt;
            najlepszaTrasa = trasa;
        }
    }

    trasaOut = najlepszaTrasa;
    return najlepszyKoszt;
}

// ---------------------------------------------------
// WIELOKROTNE URUCHOMIENIA Z RÓŻNYMI SEEDAMI
// ---------------------------------------------------
ACOResult multiRunACO(const vector<vector<double>>& D, int n, bool verbose = true) {

    // ============================================
    // KROK 1: ALGORYTM ZACHLANNY - wartość bazowa
    // ============================================
    vector<int> trasaZachlanna;
    double kosztZachlanny = algorytmZachlanny(n, D, trasaZachlanna);

    if (verbose) {
        cout << "\n=== ALGORYTM ZACHLANNY (baseline) ===" << endl;
        cout << "Koszt zachlanny (NN): " << fixed << setprecision(2) << kosztZachlanny << endl;
    }

    // Inicjalizacja globalnego wyniku wartością zachłanną BEZ 2-opt
    // ACO musi samo znaleźć lepszy wynik
    vector<int> globalNajlepszaTrasa = trasaZachlanna;
    double globalNajlepszyKoszt = kosztZachlanny;

    // ============================================
    // KROK 2: PARAMETRY ACO
    // ============================================
    ACOParams par = getAdaptiveParams(n);

    // Dla małych instancji: więcej iteracji
    if (n <= 150) {
        par.iteracje = 800;
    } else if (n <= 300) {
        par.iteracje = 500;
    }

    if (verbose) {
        cout << "\n=== PARAMETRY ADAPTACYJNE (n=" << n << ") ===\n";
        cout << "  mrowek = " << par.liczbaMrowek << endl;
        cout << "  iteracje = " << par.iteracje << endl;
        cout << "  beta = " << par.beta << endl;
        cout << "  rho = " << par.rho << endl;
        cout << "  eliteFactor = " << par.eliteFactor << endl;
        cout << "  candidateSize = " << par.candidateSize << endl;
        cout << "  opt2Passes = " << par.opt2Passes << endl;
    }

    // ============================================
    // KROK 3: WIELOKROTNE RESTARTY ACO
    // ============================================
    int maxRestarty = (n <= 150) ? 999999 : (n <= 300) ? 20 : (n <= 600) ? 5 : 3;

    int maxRestartyBezPoprawy;
    if (n <= 60) {
        maxRestartyBezPoprawy = 8;
    } else if (n <= 150) {
        maxRestartyBezPoprawy = 6;
    } else if (n <= 300) {
        maxRestartyBezPoprawy = 5;
    } else if (n <= 600) {
        maxRestartyBezPoprawy = 3;
    } else {
        maxRestartyBezPoprawy = 3;
    }

    int restartyBezPoprawy = 0;

    vector<unsigned int> seeds = {12345, 54321, 11111, 22222, 33333, 44444, 55555, 66666, 77777, 88888,
                                   99999, 13579, 24680, 36912, 48260, 59371, 60482, 71593, 82604, 93715};

    int restart = 0;
    while (restart < maxRestarty && !czasPrzekroczony()) {

        if (verbose) {
            if (n <= 150) {
                cout << "\n=== RESTART " << (restart + 1)
                     << " (pozostalo " << pozostalyCzas() << "s, bez poprawy: "
                     << restartyBezPoprawy << "/" << maxRestartyBezPoprawy << ") ===\n";
            } else {
                cout << "\n=== RESTART " << (restart + 1) << "/" << maxRestarty
                     << " (pozostalo " << pozostalyCzas() << "s, bez poprawy: "
                     << restartyBezPoprawy << "/" << maxRestartyBezPoprawy << ") ===\n";
            }
        }

        ACOResult res = runACO(D, n, par, seeds[restart % seeds.size()], verbose);

        if (res.najlepszyKoszt < globalNajlepszyKoszt - 1e-9) {
            globalNajlepszyKoszt = res.najlepszyKoszt;
            globalNajlepszaTrasa = res.najlepszaTrasa;
            restartyBezPoprawy = 0;
            if (verbose) {
                cout << "*** NOWY GLOBALNY BEST: " << fixed << setprecision(2) << globalNajlepszyKoszt << " ***\n";
            }
        } else {
            restartyBezPoprawy++;
            if (verbose) {
                cout << "Brak poprawy w tym restarcie (" << restartyBezPoprawy << "/" << maxRestartyBezPoprawy << ")\n";
            }
        }

        if (restartyBezPoprawy >= maxRestartyBezPoprawy) {
            if (verbose) {
                cout << "\n*** WCZESNE ZATRZYMANIE: brak poprawy przez " << maxRestartyBezPoprawy
                     << " restartow - prawdopodobnie osiagnieto optimum lokalne ***\n";
            }
            break;
        }

        restart++;
    }

    // ============================================
    // KROK 4: FINALNA OPTYMALIZACJA
    // ============================================
    if (verbose) cout << "\n=== FINALNA OPTYMALIZACJA ===\n";

    if (n <= 300) {
        twoOptFull(globalNajlepszaTrasa, D);
        threeOptSegment(globalNajlepszaTrasa, D, 200);
    } else {
        twoOptFast(globalNajlepszaTrasa, D, 500);
    }
    globalNajlepszyKoszt = policzKosztTrasy(globalNajlepszaTrasa, D);

    if (verbose) {
        cout << "Koszt po finalnej optymalizacji: " << fixed << setprecision(2) << globalNajlepszyKoszt << endl;
    }

    ACOResult finalRes;
    finalRes.najlepszaTrasa = globalNajlepszaTrasa;
    finalRes.najlepszyKoszt = globalNajlepszyKoszt;
    return finalRes;
}

// ---------------------------------------------------
// WCZYTYWANIE DANYCH
// ---------------------------------------------------
bool wczytajDane(const string& nazwaPliku, int& n, vector<vector<double>>& D) {
    ifstream plik(nazwaPliku);
    if (!plik.is_open()) {
        cout << "Blad otwarcia pliku: " << nazwaPliku << endl;
        return false;
    }

    plik >> n;

    string linia;
    getline(plik, linia);

    vector<pair<double, double>> punkty(n + 1);

    while (getline(plik, linia)) {
        stringstream s(linia);
        int id; double x, y;
        if (s >> id >> x >> y)
            punkty[id] = { x, y };
    }

    plik.close();

    // Macierz odległości
    D.assign(n + 1, vector<double>(n + 1, 0));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            if (i != j)
                D[i][j] = sqrt(pow(punkty[i].first - punkty[j].first, 2) +
                               pow(punkty[i].second - punkty[j].second, 2));

    return true;
}

// ---------------------------------------------------
// MAIN
// ---------------------------------------------------
int main() {
    cout << "=== ZOPTYMALIZOWANY ALGORYTM MROWKOWY (ACO) DLA TSP ===" << endl;
    cout << "Limit czasu: " << MAX_TIME_SECONDS << " sekund\n" << endl;

    cout << "Wybierz plik z danymi:" << endl;
    cout << "1. berlin52.txt" << endl;
    cout << "2. bier127.txt" << endl;
    cout << "3. tsp250.txt" << endl;
    cout << "4. tsp500.txt" << endl;
    cout << "5. tsp1000.txt" << endl;
    cout << "6. Podaj wlasna nazwe pliku" << endl;
    cout << "Wybor: ";

    int wybor;
    cin >> wybor;

    string nazwaPliku;
    switch (wybor) {
        case 1: nazwaPliku = "berlin52.txt"; break;
        case 2: nazwaPliku = "bier127.txt"; break;
        case 3: nazwaPliku = "tsp250.txt"; break;
        case 4: nazwaPliku = "tsp500.txt"; break;
        case 5: nazwaPliku = "tsp1000.txt"; break;
        case 6:
            cout << "Podaj nazwe pliku: ";
            cin >> nazwaPliku;
            break;
        default:
            cout << "Nieprawidlowy wybor!" << endl;
            return 1;
    }

    int n;
    vector<vector<double>> D;

    if (!wczytajDane(nazwaPliku, n, D)) {
        return 1;
    }

    cout << "\nWczytano " << n << " miast z pliku: " << nazwaPliku << endl;

    // START POMIARU CZASU
    globalStartTime = chrono::steady_clock::now();

    // URUCHOMIENIE ALGORYTMU
    ACOResult wynik = multiRunACO(D, n, true);

    // KONIEC POMIARU CZASU
    auto endTime = chrono::steady_clock::now();
    auto czasMs = chrono::duration_cast<chrono::milliseconds>(endTime - globalStartTime).count();

    cout << "\n========================================" << endl;
    cout << "=== WYNIK KONCOWY ===" << endl;
    cout << "========================================" << endl;
    cout << "Plik: " << nazwaPliku << endl;
    cout << "Liczba miast: " << n << endl;
    cout << "Najlepszy koszt: " << fixed << setprecision(2) << wynik.najlepszyKoszt << endl;
    cout << "Czas wykonania: " << (czasMs / 1000.0) << " s" << endl;
    cout << "========================================" << endl;

    // Opcjonalnie: wypisz trasę dla małych instancji
    if (n <= 60) {
        cout << "\nTrasa: ";
        for (int x : wynik.najlepszaTrasa) cout << x << " ";
        cout << wynik.najlepszaTrasa[0] << endl;
    }

    return 0;
}

