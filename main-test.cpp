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
#include <map>

using namespace std;

// ---------------------------------------------------
// GLOBALNE DANE DO WIZUALIZACJI
// ---------------------------------------------------
vector<pair<double, double>> globalPunkty;
string globalNazwaPliku;

// ---------------------------------------------------
// GLOBALNE ZMIENNE CZASOWE
// ---------------------------------------------------
chrono::steady_clock::time_point globalStartTime;
const int MAX_TIME_SECONDS = 30; // Krótszy czas dla testów

bool czasPrzekroczony() {
    auto teraz = chrono::steady_clock::now();
    auto sek = chrono::duration_cast<chrono::seconds>(teraz - globalStartTime).count();
    return sek >= MAX_TIME_SECONDS;
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
void twoOptFast(vector<int>& trasa, const vector<vector<double>>& D, int maxPasses = 50) {
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
    int opt2Passes;
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

    p.alfa = 1.0;
    p.Q = 100.0;
    p.tauMin = 0.001;
    p.tauMax = 20.0;

    if (n <= 100) {
        p.liczbaMrowek = min(n, 30);
        p.iteracje = 200;
        p.beta = 2.5;
        p.rho = 0.6;
        p.eliteFactor = 2;
        p.candidateSize = min(50, n - 1);
        p.opt2Passes = 50;
    }
    else if (n <= 350) {
        p.liczbaMrowek = 25;
        p.iteracje = 150;
        p.beta = 2.25;
        p.rho = 0.55;
        p.eliteFactor = 3;
        p.candidateSize = 40;
        p.opt2Passes = 30;
    }
    else {
        p.liczbaMrowek = 20;
        p.iteracje = 100;
        p.beta = 2.0;
        p.rho = 0.5;
        p.eliteFactor = 4;
        p.candidateSize = 30;
        p.opt2Passes = 20;
    }

    return p;
}

// ---------------------------------------------------
// GŁÓWNY ALGORYTM ACO (uproszczony dla testów)
// ---------------------------------------------------
ACOResult runACO(const vector<vector<double>>& D, int n, const ACOParams& par,
                 unsigned int seed, const vector<int>& trasaStartowa, double kosztStartowy) {

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

    // Lista kandydatów
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
    double tauBoost = 5.0;
    vector<vector<double>> feromony(n + 1, vector<double>(n + 1, tau0));

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

    // Atrakcyjność
    vector<vector<double>> atrakcyjnosc(n + 1, vector<double>(n + 1, 0.0));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            if (i != j) atrakcyjnosc[i][j] = 1.0 / (D[i][j] + 0.0001);

    vector<int> najlepszaTrasa = trasaStartowa;
    double najlepszyKoszt = kosztStartowy;

    int bezPoprawy = 0;
    int maxBezPoprawy = max(20, iteracje / 3);

    vector<vector<double>> tauPow(n + 1, vector<double>(n + 1));
    vector<vector<double>> etaPow(n + 1, vector<double>(n + 1));

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            etaPow[i][j] = pow(atrakcyjnosc[i][j], beta);
        }
    }

    for (int iter = 0; iter < iteracje; iter++) {

        if (czasPrzekroczony()) break;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                tauPow[i][j] = pow(feromony[i][j], alfa);
            }
        }

        vector<vector<int>> trasy(liczbaMrowek);
        vector<double> koszty(liczbaMrowek);
        bool poprawaWTejIteracji = false;

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

                for (int j : kandydaci[aktualne]) {
                    if (!odwiedzony[j]) {
                        double p = tauPow[aktualne][j] * etaPow[aktualne][j];
                        opcje.push_back({p, j});
                        suma += p;
                    }
                }

                if (suma < 1e-10) {
                    for (int j = 1; j <= n; j++) {
                        if (!odwiedzony[j]) {
                            double p = tauPow[aktualne][j] * etaPow[aktualne][j];
                            opcje.push_back({p, j});
                            suma += p;
                        }
                    }
                }

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

        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                feromony[i][j] *= (1.0 - rho);

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

        double deltaElite = eliteFactor * (Q / najlepszyKoszt);
        for (int i = 0; i < n - 1; i++) {
            int a = najlepszaTrasa[i];
            int b = najlepszaTrasa[i + 1];
            feromony[a][b] += deltaElite;
            feromony[b][a] += deltaElite;
        }
        feromony[najlepszaTrasa[n-1]][najlepszaTrasa[0]] += deltaElite;
        feromony[najlepszaTrasa[0]][najlepszaTrasa[n-1]] += deltaElite;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (feromony[i][j] < tauMin) feromony[i][j] = tauMin;
                if (feromony[i][j] > tauMax) feromony[i][j] = tauMax;
            }
        }

        if (poprawaWTejIteracji) bezPoprawy = 0;
        else bezPoprawy++;

        if (bezPoprawy >= maxBezPoprawy) break;
    }

    ACOResult res;
    res.najlepszaTrasa = najlepszaTrasa;
    res.najlepszyKoszt = najlepszyKoszt;
    return res;
}

// ---------------------------------------------------
// ALGORYTM ZACHŁANNY
// ---------------------------------------------------
double algorytmZachlanny(int n, const vector<vector<double>>& D, vector<int>& trasaOut) {
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
// GENERATOR DANYCH
// ---------------------------------------------------
void generujMiasta(int liczbaMiast, const string& nazwaPliku, unsigned int seed) {
    ofstream plikWyjsciowy(nazwaPliku);

    if (!plikWyjsciowy.is_open()) {
        cout << "Blad tworzenia pliku!" << endl;
        return;
    }

    srand(seed);

    plikWyjsciowy << liczbaMiast << endl;

    for (int i = 1; i <= liczbaMiast; i++) {
        int x = rand() % 2000 + 1;
        int y = rand() % 2500 + 1;
        plikWyjsciowy << i << " " << x << " " << y << endl;
    }

    plikWyjsciowy.close();
}

// ---------------------------------------------------
// WCZYTYWANIE DANYCH
// ---------------------------------------------------
bool wczytajDane(const string& nazwaPliku, int& n, vector<vector<double>>& D) {
    ifstream plik(nazwaPliku);
    if (!plik.is_open()) {
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

    globalPunkty = punkty;

    D.assign(n + 1, vector<double>(n + 1, 0));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            if (i != j)
                D[i][j] = sqrt(pow(punkty[i].first - punkty[j].first, 2) +
                               pow(punkty[i].second - punkty[j].second, 2));

    return true;
}

// ---------------------------------------------------
// GENEROWANIE WYKRESÓW SVG
// ---------------------------------------------------
void generujWykresSlupkowy(const string& nazwaPliku,
                           const vector<string>& etykiety,
                           const vector<double>& wartosci1,
                           const vector<double>& wartosci2,
                           const string& tytul,
                           const string& label1,
                           const string& label2,
                           const string& jednostka) {

    int width = 1000;
    int height = 600;
    int marginLeft = 80;
    int marginRight = 50;
    int marginTop = 80;
    int marginBottom = 100;

    int plotWidth = width - marginLeft - marginRight;
    int plotHeight = height - marginTop - marginBottom;

    // Znajdź maksymalną wartość
    double maxVal = 0;
    for (auto v : wartosci1) maxVal = max(maxVal, v);
    for (auto v : wartosci2) maxVal = max(maxVal, v);
    maxVal *= 1.1; // 10% margines

    ofstream f(nazwaPliku);
    if (!f.is_open()) return;

    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height << "\">\n";
    f << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Tytuł
    f << "<text x=\"" << width/2 << "\" y=\"40\" text-anchor=\"middle\" font-size=\"20\" font-weight=\"bold\" fill=\"#333\">"
      << tytul << "</text>\n";

    // Osie
    f << "<line x1=\"" << marginLeft << "\" y1=\"" << marginTop
      << "\" x2=\"" << marginLeft << "\" y2=\"" << (height - marginBottom)
      << "\" stroke=\"#333\" stroke-width=\"2\"/>\n";
    f << "<line x1=\"" << marginLeft << "\" y1=\"" << (height - marginBottom)
      << "\" x2=\"" << (width - marginRight) << "\" y2=\"" << (height - marginBottom)
      << "\" stroke=\"#333\" stroke-width=\"2\"/>\n";

    // Etykieta osi Y
    f << "<text x=\"20\" y=\"" << (marginTop + plotHeight/2)
      << "\" text-anchor=\"middle\" font-size=\"14\" fill=\"#333\" "
      << "transform=\"rotate(-90, 20, " << (marginTop + plotHeight/2) << ")\">"
      << jednostka << "</text>\n";

    int n = etykiety.size();
    double groupWidth = (double)plotWidth / n;
    double barWidth = groupWidth / 3;

    // Słupki
    for (int i = 0; i < n; i++) {
        double x = marginLeft + i * groupWidth;

        // Słupek 1 (niebieski - zachłanny)
        double h1 = (wartosci1[i] / maxVal) * plotHeight;
        double y1 = height - marginBottom - h1;
        f << "<rect x=\"" << (x + barWidth * 0.2) << "\" y=\"" << y1
          << "\" width=\"" << barWidth << "\" height=\"" << h1
          << "\" fill=\"#3366CC\" stroke=\"#1a3d7a\" stroke-width=\"1\"/>\n";

        // Słupek 2 (czerwony - mrówkowy)
        double h2 = (wartosci2[i] / maxVal) * plotHeight;
        double y2 = height - marginBottom - h2;
        f << "<rect x=\"" << (x + barWidth * 1.3) << "\" y=\"" << y2
          << "\" width=\"" << barWidth << "\" height=\"" << h2
          << "\" fill=\"#DC3912\" stroke=\"#8a2307\" stroke-width=\"1\"/>\n";

        // Etykieta osi X
        f << "<text x=\"" << (x + groupWidth/2) << "\" y=\"" << (height - marginBottom + 25)
          << "\" text-anchor=\"middle\" font-size=\"11\" fill=\"#333\">"
          << etykiety[i] << "</text>\n";
    }

    // Legenda
    int legendY = height - marginBottom + 60;
    f << "<rect x=\"" << (width/2 - 150) << "\" y=\"" << legendY
      << "\" width=\"20\" height=\"15\" fill=\"#3366CC\" stroke=\"#1a3d7a\"/>\n";
    f << "<text x=\"" << (width/2 - 125) << "\" y=\"" << (legendY + 12)
      << "\" font-size=\"13\" fill=\"#333\">" << label1 << "</text>\n";

    f << "<rect x=\"" << (width/2 + 20) << "\" y=\"" << legendY
      << "\" width=\"20\" height=\"15\" fill=\"#DC3912\" stroke=\"#8a2307\"/>\n";
    f << "<text x=\"" << (width/2 + 45) << "\" y=\"" << (legendY + 12)
      << "\" font-size=\"13\" fill=\"#333\">" << label2 << "</text>\n";

    // Linie siatki i wartości na osi Y
    for (int i = 0; i <= 5; i++) {
        double val = (maxVal / 5) * i;
        double y = height - marginBottom - (plotHeight / 5.0) * i;

        f << "<line x1=\"" << marginLeft << "\" y1=\"" << y
          << "\" x2=\"" << (width - marginRight) << "\" y2=\"" << y
          << "\" stroke=\"#ddd\" stroke-width=\"1\" stroke-dasharray=\"2,2\"/>\n";

        f << "<text x=\"" << (marginLeft - 10) << "\" y=\"" << (y + 5)
          << "\" text-anchor=\"end\" font-size=\"11\" fill=\"#666\">"
          << fixed << setprecision(0) << val << "</text>\n";
    }

    f << "</svg>\n";
    f.close();
}

void generujWykresProcentowy(const string& nazwaPliku,
                            const vector<string>& etykiety,
                            const vector<double>& procenty1,
                            const vector<double>& procenty2,
                            const string& tytul,
                            const string& label1,
                            const string& label2) {

    int width = 1000;
    int height = 600;
    int marginLeft = 80;
    int marginRight = 50;
    int marginTop = 80;
    int marginBottom = 100;

    int plotWidth = width - marginLeft - marginRight;
    int plotHeight = height - marginTop - marginBottom;

    ofstream f(nazwaPliku);
    if (!f.is_open()) return;

    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height << "\">\n";
    f << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Tytuł
    f << "<text x=\"" << width/2 << "\" y=\"40\" text-anchor=\"middle\" font-size=\"20\" font-weight=\"bold\" fill=\"#333\">"
      << tytul << "</text>\n";

    // Osie
    f << "<line x1=\"" << marginLeft << "\" y1=\"" << marginTop
      << "\" x2=\"" << marginLeft << "\" y2=\"" << (height - marginBottom)
      << "\" stroke=\"#333\" stroke-width=\"2\"/>\n";
    f << "<line x1=\"" << marginLeft << "\" y1=\"" << (height - marginBottom)
      << "\" x2=\"" << (width - marginRight) << "\" y2=\"" << (height - marginBottom)
      << "\" stroke=\"#333\" stroke-width=\"2\"/>\n";

    // Etykieta osi Y
    f << "<text x=\"20\" y=\"" << (marginTop + plotHeight/2)
      << "\" text-anchor=\"middle\" font-size=\"14\" fill=\"#333\" "
      << "transform=\"rotate(-90, 20, " << (marginTop + plotHeight/2) << ")\">Procentowa zaleznosc (%)</text>\n";

    int n = etykiety.size();
    double groupWidth = (double)plotWidth / n;
    double barWidth = groupWidth / 3;

    // Słupki (100% = maksimum)
    for (int i = 0; i < n; i++) {
        double x = marginLeft + i * groupWidth;

        // Słupek 1 (niebieski - zachłanny = 100%)
        double h1 = (procenty1[i] / 100.0) * plotHeight;
        double y1 = height - marginBottom - h1;
        f << "<rect x=\"" << (x + barWidth * 0.2) << "\" y=\"" << y1
          << "\" width=\"" << barWidth << "\" height=\"" << h1
          << "\" fill=\"#3366CC\" stroke=\"#1a3d7a\" stroke-width=\"1\"/>\n";

        // Słupek 2 (czerwony - mrówkowy)
        double h2 = (procenty2[i] / 100.0) * plotHeight;
        double y2 = height - marginBottom - h2;
        f << "<rect x=\"" << (x + barWidth * 1.3) << "\" y=\"" << y2
          << "\" width=\"" << barWidth << "\" height=\"" << h2
          << "\" fill=\"#DC3912\" stroke=\"#8a2307\" stroke-width=\"1\"/>\n";

        // Etykieta osi X
        f << "<text x=\"" << (x + groupWidth/2) << "\" y=\"" << (height - marginBottom + 25)
          << "\" text-anchor=\"middle\" font-size=\"11\" fill=\"#333\">"
          << etykiety[i] << "</text>\n";
    }

    // Legenda
    int legendY = height - marginBottom + 60;
    f << "<rect x=\"" << (width/2 - 150) << "\" y=\"" << legendY
      << "\" width=\"20\" height=\"15\" fill=\"#3366CC\" stroke=\"#1a3d7a\"/>\n";
    f << "<text x=\"" << (width/2 - 125) << "\" y=\"" << (legendY + 12)
      << "\" font-size=\"13\" fill=\"#333\">" << label1 << "</text>\n";

    f << "<rect x=\"" << (width/2 + 20) << "\" y=\"" << legendY
      << "\" width=\"20\" height=\"15\" fill=\"#DC3912\" stroke=\"#8a2307\"/>\n";
    f << "<text x=\"" << (width/2 + 45) << "\" y=\"" << (legendY + 12)
      << "\" font-size=\"13\" fill=\"#333\">" << label2 << "</text>\n";

    // Linie siatki i wartości na osi Y (80-100%)
    for (int i = 0; i <= 5; i++) {
        double val = 80.0 + (20.0 / 5.0) * i;
        double yPos = height - marginBottom - (plotHeight * ((val - 80.0) / 20.0));

        f << "<line x1=\"" << marginLeft << "\" y1=\"" << yPos
          << "\" x2=\"" << (width - marginRight) << "\" y2=\"" << yPos
          << "\" stroke=\"#ddd\" stroke-width=\"1\" stroke-dasharray=\"2,2\"/>\n";

        f << "<text x=\"" << (marginLeft - 10) << "\" y=\"" << (yPos + 5)
          << "\" text-anchor=\"end\" font-size=\"11\" fill=\"#666\">"
          << fixed << setprecision(1) << val << "</text>\n";
    }

    f << "</svg>\n";
    f.close();
}

void generujWykresBlad(const string& nazwaPliku,
                       const vector<string>& etykiety,
                       const vector<double>& bledy,
                       const string& tytul) {

    int width = 1000;
    int height = 600;
    int marginLeft = 100;
    int marginRight = 50;
    int marginTop = 80;
    int marginBottom = 100;

    int plotWidth = width - marginLeft - marginRight;
    int plotHeight = height - marginTop - marginBottom;

    // Znajdź maksymalny błąd
    double maxBlad = 0;
    for (auto b : bledy) maxBlad = max(maxBlad, b);
    maxBlad = ceil(maxBlad) + 1;

    ofstream f(nazwaPliku);
    if (!f.is_open()) return;

    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height << "\">\n";
    f << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Tytuł
    f << "<text x=\"" << width/2 << "\" y=\"40\" text-anchor=\"middle\" font-size=\"20\" font-weight=\"bold\" fill=\"#333\">"
      << tytul << "</text>\n";

    // Osie
    f << "<line x1=\"" << marginLeft << "\" y1=\"" << marginTop
      << "\" x2=\"" << marginLeft << "\" y2=\"" << (height - marginBottom)
      << "\" stroke=\"#333\" stroke-width=\"2\"/>\n";
    f << "<line x1=\"" << marginLeft << "\" y1=\"" << (height - marginBottom)
      << "\" x2=\"" << (width - marginRight) << "\" y2=\"" << (height - marginBottom)
      << "\" stroke=\"#333\" stroke-width=\"2\"/>\n";

    // Etykieta osi Y
    f << "<text x=\"20\" y=\"" << (marginTop + plotHeight/2)
      << "\" text-anchor=\"middle\" font-size=\"14\" fill=\"#333\" "
      << "transform=\"rotate(-90, 20, " << (marginTop + plotHeight/2) << ")\">Blad (%)</text>\n";

    int n = etykiety.size();
    double barWidth = (double)plotWidth / (n * 1.5);
    double spacing = barWidth * 1.5;

    // Słupki
    for (int i = 0; i < n; i++) {
        double x = marginLeft + i * spacing + spacing * 0.25;

        double h = (bledy[i] / maxBlad) * plotHeight;
        double y = height - marginBottom - h;

        f << "<rect x=\"" << x << "\" y=\"" << y
          << "\" width=\"" << barWidth << "\" height=\"" << h
          << "\" fill=\"#FF6B6B\" stroke=\"#C92A2A\" stroke-width=\"1\"/>\n";

        // Wartość na słupku
        f << "<text x=\"" << (x + barWidth/2) << "\" y=\"" << (y - 5)
          << "\" text-anchor=\"middle\" font-size=\"11\" font-weight=\"bold\" fill=\"#333\">"
          << fixed << setprecision(2) << bledy[i] << "</text>\n";

        // Etykieta osi X
        f << "<text x=\"" << (x + barWidth/2) << "\" y=\"" << (height - marginBottom + 20)
          << "\" text-anchor=\"end\" font-size=\"11\" fill=\"#333\" "
          << "transform=\"rotate(-45, " << (x + barWidth/2) << ", " << (height - marginBottom + 20) << ")\">"
          << etykiety[i] << "</text>\n";
    }

    // Linie siatki i wartości na osi Y
    int numLines = min(10, (int)maxBlad + 1);
    for (int i = 0; i <= numLines; i++) {
        double val = (maxBlad / numLines) * i;
        double y = height - marginBottom - (plotHeight / (double)numLines) * i;

        f << "<line x1=\"" << marginLeft << "\" y1=\"" << y
          << "\" x2=\"" << (width - marginRight) << "\" y2=\"" << y
          << "\" stroke=\"#ddd\" stroke-width=\"1\" stroke-dasharray=\"2,2\"/>\n";

        f << "<text x=\"" << (marginLeft - 10) << "\" y=\"" << (y + 5)
          << "\" text-anchor=\"end\" font-size=\"11\" fill=\"#666\">"
          << fixed << setprecision(1) << val << "</text>\n";
    }

    f << "</svg>\n";
    f.close();
}

// ---------------------------------------------------
// TEST 1: Porównanie ACO z algorytmem zachłannym
// ---------------------------------------------------
void test1_PorownanieZZachlannymLosowe() {
    cout << "\n========================================" << endl;
    cout << "TEST 1: Porownanie ACO z algorytmem zachlannym" << endl;
    cout << "15 losowych instancji (40-200 miast)" << endl;
    cout << "========================================\n" << endl;

    vector<int> rozmiary;
    // 15 punktów pomiarowych z zakresu 40-200, w miarę równymi odstępami
    for (int i = 0; i < 15; i++) {
        rozmiary.push_back(40 + i * 11); // 40, 51, 62, 73, 84, 95, 106, 117, 128, 139, 150, 161, 172, 183, 194
    }

    vector<string> etykiety;
    vector<double> wyniki_zachlanny;
    vector<double> wyniki_aco;
    vector<double> procenty_zachlanny;
    vector<double> procenty_aco;

    for (int i = 0; i < rozmiary.size(); i++) {
        int n = rozmiary[i];
        string nazwaPliku = "test_losowy_" + to_string(n) + ".txt";

        cout << "Test " << (i+1) << "/15: " << n << " miast..." << flush;

        // Generuj dane
        generujMiasta(n, nazwaPliku, 12345 + i * 100);

        // Wczytaj dane
        vector<vector<double>> D;
        int nCzyt;
        if (!wczytajDane(nazwaPliku, nCzyt, D)) {
            cout << " BLAD wczytywania!" << endl;
            continue;
        }

        globalStartTime = chrono::steady_clock::now();

        // Algorytm zachłanny
        vector<int> trasaZachlanna;
        double kosztZachlanny = algorytmZachlanny(n, D, trasaZachlanna);

        // ACO
        ACOParams par = getAdaptiveParams(n);
        ACOResult resACO = runACO(D, n, par, 12345, trasaZachlanna, kosztZachlanny);

        cout << " Zachlanny: " << fixed << setprecision(0) << kosztZachlanny
             << ", ACO: " << kosztZachlanny << " -> " << resACO.najlepszyKoszt << endl;

        etykiety.push_back(to_string(n));
        wyniki_zachlanny.push_back(kosztZachlanny);
        wyniki_aco.push_back(resACO.najlepszyKoszt);
        procenty_zachlanny.push_back(100.0);
        procenty_aco.push_back((resACO.najlepszyKoszt / kosztZachlanny) * 100.0);
    }

    // Generuj wykresy
    cout << "\nGenerowanie wykresow..." << endl;

    generujWykresSlupkowy("Wykres1Dystans.svg", etykiety, wyniki_zachlanny, wyniki_aco,
                          "Algorytm zachlanny a algorytm mrowkowy (Dystans)",
                          "zachlanny", "mrowkowy", "Dystans");
    cout << "Zapisano: Wykres1Dystans.svg" << endl;

    generujWykresProcentowy("Wykres1Procent.svg", etykiety, procenty_zachlanny, procenty_aco,
                           "Algorytm zachlanny a mrowkowy (w procentach)",
                           "zachlanny", "mrowkowy");
    cout << "Zapisano: Wykres1Procent.svg" << endl;
}

// ---------------------------------------------------
// TEST 2: Porównanie z wartościami optymalnymi (benchmarki)
// ---------------------------------------------------
void test2_PorownanieZOptimum() {
    cout << "\n========================================" << endl;
    cout << "TEST 2: Porownanie z wartosciami optymalnymi" << endl;
    cout << "10 instancji benchmarkowych" << endl;
    cout << "========================================\n" << endl;

    map<string, double> optymalne;
    optymalne["berlin52"] = 7542;
    optymalne["eil51"] = 426;
    optymalne["eil76"] = 538;
    optymalne["st70"] = 675;
    optymalne["pr76"] = 108159;
    optymalne["kroA100"] = 21282;
    optymalne["kroB100"] = 22141;
    optymalne["kroC100"] = 20749;
    optymalne["kroD100"] = 21294;
    optymalne["lin105"] = 14379;

    vector<string> instancje = {"berlin52", "eil51", "eil76", "st70", "pr76",
                                "kroA100", "kroB100", "kroC100", "kroD100", "lin105"};

    vector<string> etykiety;
    vector<double> procenty_aco;
    vector<double> procenty_opt;
    vector<double> bledy;

    for (int i = 0; i < instancje.size(); i++) {
        string nazwa = instancje[i];
        string nazwaPliku = nazwa + ".txt";

        cout << "Test " << (i+1) << "/10: " << nazwa << "..." << flush;

        // Wczytaj dane
        vector<vector<double>> D;
        int n;
        if (!wczytajDane(nazwaPliku, n, D)) {
            cout << " BLAD - brak pliku!" << endl;
            continue;
        }

        globalStartTime = chrono::steady_clock::now();

        // Algorytm zachłanny
        vector<int> trasaZachlanna;
        double kosztZachlanny = algorytmZachlanny(n, D, trasaZachlanna);

        // ACO
        ACOParams par = getAdaptiveParams(n);
        ACOResult resACO = runACO(D, n, par, 54321, trasaZachlanna, kosztZachlanny);

        double opt = optymalne[nazwa];
        double blad = ((resACO.najlepszyKoszt - opt) / opt) * 100.0;

        cout << " Optimum: " << fixed << setprecision(0) << opt
             << ", ACO: " << resACO.najlepszyKoszt
             << " (blad: " << setprecision(2) << blad << "%)" << endl;

        etykiety.push_back(nazwa);
        procenty_aco.push_back((resACO.najlepszyKoszt / opt) * 100.0);
        procenty_opt.push_back(100.0);
        bledy.push_back(blad);
    }

    // Generuj wykresy
    cout << "\nGenerowanie wykresow..." << endl;

    generujWykresProcentowy("Wykres2Procent.svg", etykiety, procenty_aco, procenty_opt,
                           "Algorytm mrowkowy a optimum (w procentach)",
                           "mrowkowy", "optimum");
    cout << "Zapisano: Wykres2Procent.svg" << endl;

    generujWykresBlad("Wykres2Blad.svg", etykiety, bledy,
                     "Blad wzgledny algorytmu wzgledem wartosci optymalnej");
    cout << "Zapisano: Wykres2Blad.svg" << endl;
}

// ---------------------------------------------------
// MAIN
// ---------------------------------------------------
int main() {
    cout << "=== TESTY AUTOMATYCZNE ACO ===" << endl;
    cout << "Program przeprowadzi 2 serie testow i wygeneruje wykresy\n" << endl;

    // Test 1: Porównanie z algorytmem zachłannym (losowe instancje)
    test1_PorownanieZZachlannymLosowe();

    // Test 2: Porównanie z wartościami optymalnymi (benchmarki)
    test2_PorownanieZOptimum();

    cout << "\n========================================" << endl;
    cout << "=== TESTY ZAKONCZONE ===" << endl;
    cout << "========================================" << endl;
    cout << "\nWygenerowane wykresy:" << endl;
    cout << "1. Wykres1Dystans.svg - Porownanie dystansow (zachlanny vs ACO)" << endl;
    cout << "2. Wykres1Procent.svg - Porownanie procentowe (zachlanny vs ACO)" << endl;
    cout << "3. Wykres2Procent.svg - Porownanie z optimum (procentowo)" << endl;
    cout << "4. Wykres2Blad.svg - Blad wzgledny wzgledem optimum" << endl;
    cout << "\nPliki SVG mozna otworzyc w przegladarce internetowej." << endl;

    return 0;
}

