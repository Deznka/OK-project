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
// GLOBALNE DANE DO WIZUALIZACJI
// ---------------------------------------------------
vector<pair<double, double>> globalPunkty;
string globalNazwaPliku;

// ---------------------------------------------------
// FUNKCJE WIZUALIZACJI - GENEROWANIE SVG
// ---------------------------------------------------
void generujSVG(const string& nazwaPliku,
                const vector<pair<double, double>>& punkty,
                const vector<int>& trasa,
                const string& tytul,
                bool rysujTrase = true) {

    int n = (int)punkty.size() - 1; // punkty[0] nie używane
    if (n <= 0) return;

    // Znajdź zakres współrzędnych
    double minX = numeric_limits<double>::infinity();
    double maxX = -numeric_limits<double>::infinity();
    double minY = numeric_limits<double>::infinity();
    double maxY = -numeric_limits<double>::infinity();

    for (int i = 1; i <= n; i++) {
        minX = min(minX, punkty[i].first);
        maxX = max(maxX, punkty[i].first);
        minY = min(minY, punkty[i].second);
        maxY = max(maxY, punkty[i].second);
    }

    // Margines
    double marginX = (maxX - minX) * 0.1;
    double marginY = (maxY - minY) * 0.1;
    minX -= marginX; maxX += marginX;
    minY -= marginY; maxY += marginY;

    // Rozmiar SVG
    int svgWidth = 800;
    int svgHeight = 700;
    int margin = 50;

    // Funkcje skalowania
    auto scaleX = [&](double x) -> double {
        return margin + (x - minX) / (maxX - minX) * (svgWidth - 2 * margin);
    };
    auto scaleY = [&](double y) -> double {
        // Odwrócenie Y bo w SVG Y rośnie w dół
        return svgHeight - margin - (y - minY) / (maxY - minY) * (svgHeight - 2 * margin);
    };

    ofstream plik(nazwaPliku);
    if (!plik.is_open()) {
        cout << "Nie mozna utworzyc pliku SVG: " << nazwaPliku << endl;
        return;
    }

    // Nagłówek SVG
    plik << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    plik << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << svgWidth << "\" height=\"" << svgHeight + 50 << "\">\n";

    // Tło
    plik << "  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Tytuł
    plik << "  <text x=\"" << svgWidth/2 << "\" y=\"30\" text-anchor=\"middle\" font-size=\"18\" font-weight=\"bold\" fill=\"#333\">"
         << tytul << "</text>\n";

    // Osie współrzędnych (uproszczone)
    plik << "  <line x1=\"" << margin << "\" y1=\"" << svgHeight - margin
         << "\" x2=\"" << svgWidth - margin << "\" y2=\"" << svgHeight - margin
         << "\" stroke=\"#ccc\" stroke-width=\"1\"/>\n";
    plik << "  <line x1=\"" << margin << "\" y1=\"" << margin
         << "\" x2=\"" << margin << "\" y2=\"" << svgHeight - margin
         << "\" stroke=\"#ccc\" stroke-width=\"1\"/>\n";

    // Etykiety osi
    plik << "  <text x=\"" << svgWidth/2 << "\" y=\"" << svgHeight - 10
         << "\" text-anchor=\"middle\" font-size=\"12\" fill=\"#666\">X</text>\n";
    plik << "  <text x=\"15\" y=\"" << svgHeight/2
         << "\" text-anchor=\"middle\" font-size=\"12\" fill=\"#666\" transform=\"rotate(-90, 15, " << svgHeight/2 << ")\">Y</text>\n";

    // Rysuj trasę (linie między miastami)
    if (rysujTrase && !trasa.empty()) {
        plik << "  <!-- Trasa -->\n";
        plik << "  <path d=\"M";

        for (size_t i = 0; i < trasa.size(); i++) {
            int miasto = trasa[i];
            double sx = scaleX(punkty[miasto].first);
            double sy = scaleY(punkty[miasto].second);

            if (i == 0) {
                plik << fixed << setprecision(2) << sx << "," << sy;
            } else {
                plik << " L" << fixed << setprecision(2) << sx << "," << sy;
            }
        }

        // Zamknięcie cyklu - powrót do pierwszego miasta
        if (!trasa.empty()) {
            double sx = scaleX(punkty[trasa[0]].first);
            double sy = scaleY(punkty[trasa[0]].second);
            plik << " L" << fixed << setprecision(2) << sx << "," << sy;
        }

        plik << "\" fill=\"none\" stroke=\"#2196F3\" stroke-width=\"2\" stroke-opacity=\"0.7\"/>\n";
    }

    // Rysuj punkty (miasta)
    plik << "  <!-- Miasta -->\n";
    double r = (n <= 100) ? 5 : (n <= 300) ? 4 : 3; // Rozmiar punktu zależny od liczby miast

    for (int i = 1; i <= n; i++) {
        double sx = scaleX(punkty[i].first);
        double sy = scaleY(punkty[i].second);

        // Punkt
        plik << "  <circle cx=\"" << fixed << setprecision(2) << sx
             << "\" cy=\"" << sy << "\" r=\"" << r << "\" fill=\"#E53935\" stroke=\"#B71C1C\" stroke-width=\"1\"/>\n";

        // Etykieta (tylko dla małych instancji)
        if (n <= 60) {
            plik << "  <text x=\"" << (sx + r + 2) << "\" y=\"" << (sy - r - 2)
                 << "\" font-size=\"9\" fill=\"#333\">" << i << "</text>\n";
        }
    }

    // Oznaczenie miasta startowego (jeśli jest trasa)
    if (rysujTrase && !trasa.empty()) {
        int start = trasa[0];
        double sx = scaleX(punkty[start].first);
        double sy = scaleY(punkty[start].second);
        plik << "  <circle cx=\"" << fixed << setprecision(2) << sx
             << "\" cy=\"" << sy << "\" r=\"" << (r + 3) << "\" fill=\"none\" stroke=\"#4CAF50\" stroke-width=\"3\"/>\n";
    }

    // Legenda
    int legendaY = svgHeight + 20;
    plik << "  <circle cx=\"" << margin << "\" cy=\"" << legendaY << "\" r=\"5\" fill=\"#E53935\"/>\n";
    plik << "  <text x=\"" << (margin + 10) << "\" y=\"" << (legendaY + 4) << "\" font-size=\"11\" fill=\"#333\">Miasta (n=" << n << ")</text>\n";

    if (rysujTrase && !trasa.empty()) {
        plik << "  <line x1=\"" << (margin + 150) << "\" y1=\"" << legendaY
             << "\" x2=\"" << (margin + 180) << "\" y2=\"" << legendaY
             << "\" stroke=\"#2196F3\" stroke-width=\"2\"/>\n";
        plik << "  <text x=\"" << (margin + 185) << "\" y=\"" << (legendaY + 4) << "\" font-size=\"11\" fill=\"#333\">Trasa</text>\n";

        plik << "  <circle cx=\"" << (margin + 280) << "\" cy=\"" << legendaY << "\" r=\"8\" fill=\"none\" stroke=\"#4CAF50\" stroke-width=\"2\"/>\n";
        plik << "  <text x=\"" << (margin + 293) << "\" y=\"" << (legendaY + 4) << "\" font-size=\"11\" fill=\"#333\">Start</text>\n";
    }

    plik << "</svg>\n";
    plik.close();

    cout << "Zapisano graf do pliku: " << nazwaPliku << endl;
}

// Przeciążona wersja bez trasy (tylko punkty)
void generujSVGBezTrasy(const string& nazwaPliku,
                        const vector<pair<double, double>>& punkty,
                        const string& tytul) {
    vector<int> pustaTrasa;
    generujSVG(nazwaPliku, punkty, pustaTrasa, tytul, false);
}

// ---------------------------------------------------
// GENERATOR DANYCH
// ---------------------------------------------------
void generujMiasta(int liczbaMiast, const string& nazwaPliku) {
    ofstream plikWyjsciowy(nazwaPliku);

    if (!plikWyjsciowy.is_open()) {
        cout << "Blad tworzenia pliku!" << endl;
        return;
    }

    srand(time(0));

    plikWyjsciowy << liczbaMiast << endl;

    for (int i = 1; i <= liczbaMiast; i++) {
        int x = rand() % 2000 + 1;
        int y = rand() % 2500 + 1;
        plikWyjsciowy << i << " " << x << " " << y << endl;
    }

    plikWyjsciowy.close();
    cout << "Wygenerowano " << liczbaMiast << " miast i zapisano do pliku " << nazwaPliku << endl;
}

// ---------------------------------------------------
// GLOBALNE ZMIENNE CZASOWE
// ---------------------------------------------------
chrono::steady_clock::time_point globalStartTime;
const int MAX_TIME_SECONDS = 175;

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
    bool bylRestartFeromonow;  // czy podczas działania nastąpił restart feromonów
    vector<int> trasaPrzyRestarcie;  // trasa w momencie restartu feromonów (stagnacja)
    double kosztPrzyRestarcie;       // koszt w momencie restartu feromonów
    vector<int> trasaPo100Iteracjach;  // trasa po pierwszych 100 iteracjach
    double kosztPo100Iteracjach;       // koszt po pierwszych 100 iteracjach
    bool zapisanoPo100;                // czy zapisano trasę po 100 iteracjach
};

// ---------------------------------------------------
// PARAMETRY ADAPTACYJNE ZALEŻNE OD ROZMIARU PROBLEMU
// ---------------------------------------------------
ACOParams getAdaptiveParams(int n) {
    ACOParams p;

    // Domyślne, "bezpieczne" wartości bazowe
    p.alfa = 1.0;           // Waga feromonów - standardowo 1.0
    p.Q = 100.0;            // Stała feromonowa
    p.tauMin = 0.001;       // Dolny limit feromonu (zapobiega stagnacji w 0)
    p.tauMax = 20.0;        // Górny limit feromonu

    // Parametry adaptacyjne zależne od wielkości problemu (N)
    // Podział na koszyki wielkości, uniwersalny dla różnych instancji

    if (n <= 100) {
        // MAŁE INSTANCJE (do 100 miast)
        // Możemy pozwolić sobie na więcej mrowek i agresywniejszą eksplorację
        p.liczbaMrowek = min(n, 50);    // Sporo mrówek dla dokładnego przeszukania
        p.iteracje = 800;               // Dużo iteracji, bo są szybkie
        p.beta = 2.5;                   // Umiarkowany wpływ heurystyki (wyrównane szanse)
        p.rho = 0.6;                    // Powolne parowanie = stabilne uczenie
        p.eliteFactor = 2;              // Lekkie wzmocnienie najlepszego
        p.candidateSize = min(50, n - 1); // Szerokie przeszukiwanie sąsiedztwa
        p.opt2Passes = 100;             // Intensywne poprawianie 2-opt
    }
    else if (n <= 350) {
        // ŚREDNIE INSTANCJE (101 - 350 miast)
        // Balans między szybkością a jakością
        p.liczbaMrowek = 40;            // Stała liczba mrówek wystarczy
        p.iteracje = 500;               // Mniej iteracji
        p.beta = 2.25;                   // Większy nacisk na odległość (heurystykę)
        p.rho = 0.55;
        p.eliteFactor = 3;              // Mocniejsze promowanie lidera
        p.candidateSize = 40;           // Ograniczenie kandydatów
        p.opt2Passes = 60;              // Mniej optymalizacji lokalnej
    }
    else if (n <= 700) {
        // DUŻE INSTANCJE (351 - 700 miast)
        // Optymalizacja pod kątem wydajności
        p.liczbaMrowek = 30;
        p.iteracje = 350;
        p.beta = 2;                   // Silna heurystyka (ważna dla dużych grafów)
        p.rho = 0.4;
        p.eliteFactor = 4;
        p.candidateSize = 30;           // Węższe sąsiedztwo
        p.opt2Passes = 40;              // Jeszcze mniej 2-opt
    }
    else {
        // BARDZO DUŻE INSTANCJE (> 700 miast)
        // Priorytet: zmieścić się w czasie
        p.liczbaMrowek = 20;
        p.iteracje = 200;               // Mało iteracji
        p.beta = 3;                   // Bardzo silna heurystyka (greedy-like)
        p.rho = 0.5;
        p.eliteFactor = 5;              // Bardzo silne wzmocnienie najlepszego
        p.candidateSize = 25;           // Bardzo wąskie sąsiedztwo
        p.opt2Passes = 20;              // Minimalny local search
    }

    return p;
}

// ---------------------------------------------------
// GŁÓWNY ALGORYTM ACO
// ---------------------------------------------------
ACOResult runACO(const vector<vector<double>>& D, int n, const ACOParams& par,
                 unsigned int seed, const vector<int>& trasaStartowa, double kosztStartowy,
                 bool verbose = false) {

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

    // Śledzenie restartu feromonów (stagnacja)
    bool bylRestartFeromonow = false;
    vector<int> trasaPrzyRestarcie;
    double kosztPrzyRestarcie = 0.0;

    // Śledzenie trasy po 100 iteracjach
    vector<int> trasaPo100Iteracjach;
    double kosztPo100Iteracjach = 0.0;
    bool zapisanoPo100 = false;

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

    // Inicjalizacja feromonów - START OD TRASY ZACHŁANNEJ
    double tau0 = 1.0;
    double tauBoost = 5.0;
    vector<vector<double>> feromony(n + 1, vector<double>(n + 1, tau0));

    // Wzmocnij feromony na krawędziach trasy startowej (zachłannej)
    if (!trasaStartowa.empty()) {
        for (int i = 0; i < n - 1; i++) {
            int a = trasaStartowa[i];
            int b = trasaStartowa[i + 1];
            feromony[a][b] = tauBoost;
            feromony[b][a] = tauBoost;
        }
        feromony[trasaStartowa[n-1]][trasaStartowa[0]] = tauBoost;
        feromony[trasaStartowa[0]][trasaStartowa[n-1]] = tauBoost;
    }

    // Atrakcyjność (1/odległość)
    vector<vector<double>> atrakcyjnosc(n + 1, vector<double>(n + 1, 0.0));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            if (i != j) atrakcyjnosc[i][j] = 1.0 / (D[i][j] + 0.0001);

    // STARTUJEMY OD TRASY ZACHŁANNEJ
    vector<int> najlepszaTrasa = trasaStartowa;
    double najlepszyKoszt = kosztStartowy;

    int bezPoprawy = 0;
    int maxBezPoprawy = max(30, iteracje / 2);

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

        // === ZAPISZ TRASĘ PO 100 ITERACJACH ===
        if (iter == 100 && !zapisanoPo100) {
            trasaPo100Iteracjach = najlepszaTrasa;
            kosztPo100Iteracjach = najlepszyKoszt;
            zapisanoPo100 = true;
            if (verbose) {
                cout << ">>> Zapisano stan po 100 iteracjach: koszt = " << fixed << setprecision(2) << kosztPo100Iteracjach << endl;
            }
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

        if (verbose && (iter % 50 == 0 || poprawaWTejIteracji)) {
            cout << "Iter " << iter + 1 << " | Best: " << fixed << setprecision(2) << najlepszyKoszt << endl;
        }

        // Warunek stopu - brak poprawy
        if (poprawaWTejIteracji) bezPoprawy = 0;
        else bezPoprawy++;

        if (bezPoprawy >= maxBezPoprawy) {
            if (verbose) cout << "Stop: brak poprawy przez " << bezPoprawy << " iteracji.\n";
            break;
        }

        // MOCNIEJSZY RESTART - perturbacja zamiast pełnego resetu
        if (bezPoprawy > 0 && bezPoprawy % (maxBezPoprawy / 3) == 0) {
            if (verbose) cout << "Perturbacja feromonow (stagnacja).\n";

            // Zapisz trasę przy pierwszym restarcie feromonów (stagnacja)
            if (!bylRestartFeromonow) {
                bylRestartFeromonow = true;
                trasaPrzyRestarcie = najlepszaTrasa;
                kosztPrzyRestarcie = najlepszyKoszt;
            }

            // Zamiast pełnego resetu - dodaj szum i wzmocnij najlepszą trasę
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= n; j++) {
                    double szum = (rand() / (double)RAND_MAX) * tau0 * 0.5;
                    feromony[i][j] = feromony[i][j] * 0.5 + tau0 * 0.3 + szum;
                }
            }

            // Mocno wzmocnij najlepszą trasę po perturbacji
            for (int i = 0; i < n - 1; i++) {
                int a = najlepszaTrasa[i];
                int b = najlepszaTrasa[i + 1];
                feromony[a][b] += tauBoost;
                feromony[b][a] += tauBoost;
            }
            feromony[najlepszaTrasa[n-1]][najlepszaTrasa[0]] += tauBoost;
            feromony[najlepszaTrasa[0]][najlepszaTrasa[n-1]] += tauBoost;
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
    res.bylRestartFeromonow = bylRestartFeromonow;
    res.trasaPrzyRestarcie = trasaPrzyRestarcie;
    res.kosztPrzyRestarcie = kosztPrzyRestarcie;
    res.trasaPo100Iteracjach = trasaPo100Iteracjach;
    res.kosztPo100Iteracjach = kosztPo100Iteracjach;
    res.zapisanoPo100 = zapisanoPo100;
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

    // === WIZUALIZACJA 2: Graf po algorytmie zachłannym ===
    if (!globalPunkty.empty()) {
        stringstream ss;
        ss << "Algorytm Zachlanny (Nearest Neighbor) - Koszt: " << fixed << setprecision(2) << kosztZachlanny;
        generujSVG("graf_2_zachlanny.svg", globalPunkty, trasaZachlanna, ss.str(), true);
    }

    // Optymalizacja trasy zachłannej przez 2-opt jako punkt startowy
    vector<int> trasaZachlannaOpt = trasaZachlanna;
    twoOptFull(trasaZachlannaOpt, D);
    double kosztZachlannyOpt = policzKosztTrasy(trasaZachlannaOpt, D);

    // Inicjalizacja globalnego wyniku - startujemy od zoptymalizowanej trasy zachłannej
    vector<int> globalNajlepszaTrasa = trasaZachlannaOpt;
    double globalNajlepszyKoszt = kosztZachlannyOpt;

    // Zapamiętaj oryginalny koszt zachłanny (bez 2-opt) do porównania
    double oryginalnyKosztZachlanny = kosztZachlanny;

    // ============================================
    // KROK 2: PARAMETRY ACO
    // ============================================
    ACOParams par = getAdaptiveParams(n);

    // Zwiększ liczbę iteracji dla lepszej eksploracji
    par.iteracje = max(par.iteracje, 600);

    if (verbose) {
        cout << "\n=== PARAMETRY ACO (uniwersalne) (n=" << n << ") ===\n";
        cout << "  alfa = " << par.alfa << " (waga feromonow)" << endl;
        cout << "  beta = " << par.beta << " (waga heurystyki)" << endl;
        cout << "  rho = " << par.rho << " (parowanie feromonow)" << endl;
        cout << "  Q = " << par.Q << " (stala aktualizacji)" << endl;
        cout << "  mrowek = " << par.liczbaMrowek << endl;
        cout << "  iteracje = " << par.iteracje << endl;
        cout << "  candidateSize = " << par.candidateSize << endl;
    }

    // ============================================
    // KROK 3: WIELOKROTNE RESTARTY ACO
    // ============================================
    int maxRestartyBezPoprawy = 15;
    bool znalezionoLepszeNizZachlanny = false;
    int restartyBezPoprawy = 0;

    vector<unsigned int> seeds = {12345, 54321, 11111, 22222, 33333, 44444, 55555, 66666, 77777, 88888,
                                   99999, 13579, 24680, 36912, 48260, 59371, 60482, 71593, 82604, 93715,
                                   10101, 20202, 30303, 40404, 50505, 60606, 70707, 80808, 90909, 12121};

    // Zmienne do śledzenia wizualizacji ACO
    // Graf 3: po 100 iteracjach
    // Graf 4: przy pierwszej stagnacji (perturbacji feromonów)
    bool zapisanoACO1 = false;  // Graf 3 - po 100 iteracjach
    bool zapisanoACO2 = false;  // Graf 4 - przy stagnacji

    int restart = 0;
    while (!czasPrzekroczony()) {

        if (verbose) {
            cout << "\n=== RESTART " << (restart + 1)
                 << " (pozostalo " << pozostalyCzas() << "s, bez poprawy: "
                 << restartyBezPoprawy << "/" << maxRestartyBezPoprawy << ") ===\n";
        }

        // Zmieniaj parametry co kilka restartów dla większej różnorodności
        ACOParams parRun = par;
        if (restart % 3 == 1) {
            parRun.beta *= 1.5;  // Więcej heurystyki
        } else if (restart % 3 == 2) {
            parRun.beta *= 0.7;  // Mniej heurystyki, więcej eksploracji
            parRun.rho *= 1.5;   // Szybsze parowanie
        }

        ACOResult res = runACO(D, n, parRun, seeds[restart % seeds.size()], trasaZachlanna, kosztZachlanny, verbose);

        // === WIZUALIZACJA 3: Po 100 iteracjach ===
        if (!zapisanoACO1 && !globalPunkty.empty() && res.zapisanoPo100 && !res.trasaPo100Iteracjach.empty()) {
            stringstream ss;
            ss << "ACO - Po 100 iteracjach - Koszt: " << fixed << setprecision(2) << res.kosztPo100Iteracjach;
            generujSVG("graf_3_aco_etap1.svg", globalPunkty, res.trasaPo100Iteracjach, ss.str(), true);
            zapisanoACO1 = true;
            if (verbose) {
                cout << ">>> Zapisano graf etapu 1 (po 100 iteracjach)\n";
            }
        }

        // === WIZUALIZACJA 4: Przy pierwszej stagnacji (perturbacji feromonów) ===
        if (!zapisanoACO2 && !globalPunkty.empty() && res.bylRestartFeromonow && !res.trasaPrzyRestarcie.empty()) {
            stringstream ss;
            ss << "ACO - Pierwsza stagnacja (perturbacja) - Koszt: " << fixed << setprecision(2) << res.kosztPrzyRestarcie;
            generujSVG("graf_4_aco_etap2.svg", globalPunkty, res.trasaPrzyRestarcie, ss.str(), true);
            zapisanoACO2 = true;
            if (verbose) {
                cout << ">>> Zapisano graf etapu 2 (pierwsza stagnacja)\n";
            }
        }

        // Aktualizacja globalnego wyniku
        if (res.najlepszyKoszt < globalNajlepszyKoszt - 1e-9) {
            globalNajlepszyKoszt = res.najlepszyKoszt;
            globalNajlepszaTrasa = res.najlepszaTrasa;
            restartyBezPoprawy = 0;

            if (globalNajlepszyKoszt < oryginalnyKosztZachlanny - 1e-9) {
                znalezionoLepszeNizZachlanny = true;
            }

            if (verbose) {
                cout << "*** NOWY NAJLEPSZY: " << fixed << setprecision(2) << globalNajlepszyKoszt;
                cout << " ***\n";
            }
        } else {
            restartyBezPoprawy++;
            if (verbose) {
                cout << "Brak poprawy w tym restarcie (" << restartyBezPoprawy << "/" << maxRestartyBezPoprawy << ")\n";
            }
        }

        // Warunek zakończenia
        if (znalezionoLepszeNizZachlanny && restartyBezPoprawy >= maxRestartyBezPoprawy) {
            if (verbose) {
                cout << "\n*** ZAKONCZENIE: znaleziono lepsze niz zachlanny i brak poprawy przez "
                     << maxRestartyBezPoprawy << " restartow ***\n";
            }
            break;
        }

        if (!znalezionoLepszeNizZachlanny && restartyBezPoprawy >= maxRestartyBezPoprawy) {
            if (verbose) {
                cout << ">>> Nie znaleziono lepszego niz zachłanny - kontynuuje szukanie...\n";
            }
            restartyBezPoprawy = maxRestartyBezPoprawy / 2;
        }

        restart++;
    }

    // Jeśli nie było grafu etapu 1, zapisz aktualny stan
    if (!zapisanoACO1 && !globalPunkty.empty()) {
        stringstream ss;
        ss << "ACO - Stan poczatkowy - Koszt: " << fixed << setprecision(2) << globalNajlepszyKoszt;
        generujSVG("graf_3_aco_etap1.svg", globalPunkty, globalNajlepszaTrasa, ss.str(), true);
    }

    // Jeśli nie było grafu etapu 2, zapisz stan przed optymalizacją finalną
    if (!zapisanoACO2 && !globalPunkty.empty()) {
        stringstream ss;
        ss << "ACO - Przed optymalizacja finalna - Koszt: " << fixed << setprecision(2) << globalNajlepszyKoszt;
        generujSVG("graf_4_aco_etap2.svg", globalPunkty, globalNajlepszaTrasa, ss.str(), true);
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
        cout << "Koszt po finalnej optymalizacji: " << fixed << setprecision(2) << globalNajlepszyKoszt;
        if (globalNajlepszyKoszt < oryginalnyKosztZachlanny - 1e-9) {
            cout << " >>> SUKCES: Wynik lepszy niz algorytm zachlanny! ("
                 << fixed << setprecision(2) << (oryginalnyKosztZachlanny - globalNajlepszyKoszt)
                 << " lepiej)";
        }
        cout << endl;
    }

    // === WIZUALIZACJA 5: Finalna trasa ===
    if (!globalPunkty.empty()) {
        stringstream ss;
        ss << "WYNIK KONCOWY - Koszt: " << fixed << setprecision(2) << globalNajlepszyKoszt;
        generujSVG("graf_5_finalny.svg", globalPunkty, globalNajlepszaTrasa, ss.str(), true);
    }

    ACOResult finalRes;
    finalRes.najlepszaTrasa = globalNajlepszaTrasa;
    finalRes.najlepszyKoszt = globalNajlepszyKoszt;
    finalRes.bylRestartFeromonow = false;
    finalRes.zapisanoPo100 = false;
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

    // Zapisz punkty globalnie dla wizualizacji
    globalPunkty = punkty;

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
    cout << "7. Wygeneruj losowy graf (generator.txt)" << endl;
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
        case 7: {
            int liczba;
            cout << "Podaj liczbe miast: ";
            cin >> liczba;
            generujMiasta(liczba, "generator.txt");
            nazwaPliku = "generator.txt";
            break;
        }
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

    // === WIZUALIZACJA 1: Graf początkowy (tylko punkty, bez połączeń) ===
    if (!globalPunkty.empty()) {
        stringstream ss;
        ss << "Instancja poczatkowa - " << n << " miast (bez polaczen)";
        generujSVGBezTrasy("graf_1_poczatkowy.svg", globalPunkty, ss.str());
    }

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

    // Podsumowanie wygenerowanych grafów
    cout << "\n=== WYGENEROWANE GRAFY SVG ===" << endl;
    cout << "1. graf_1_poczatkowy.svg   - Instancja poczatkowa (tylko miasta)" << endl;
    cout << "2. graf_2_zachlanny.svg    - Trasa po algorytmie zachlannym" << endl;
    cout << "3. graf_3_aco_etap1.svg    - ACO: pierwsza znaczaca poprawa" << endl;
    cout << "4. graf_4_aco_etap2.svg    - ACO: etap posredni" << endl;
    cout << "5. graf_5_finalny.svg      - Wynik koncowy" << endl;
    cout << "\nPliki SVG mozna otworzyc w przegladarce internetowej." << endl;

    if (n <= 60) {
        cout << "\nTrasa: ";
        for (int x : wynik.najlepszaTrasa) cout << x << " ";
        cout << wynik.najlepszaTrasa[0] << endl;
    }

    return 0;
}

