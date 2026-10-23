#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <chrono>
#include <climits>
using namespace std;

ifstream fin("input.txt");
ofstream fout("output.txt");

const double INF = 1e18;
const long long BIG = LLONG_MAX / 4;
const size_t MAXSEQ = 5000000;      // limita de siguranta pentru lungimea secventei
const double TIME_PER_CASE = 5.0;   // secunde de cautare locala per caz

typedef vector<pair<int, long long>> Mols; // (id molecula, cantitate)

struct Reaction {
    string text;              // linia originala, pentru output
    Mols in, out;
    vector<long long> outNet; // productia neta pentru fiecare iesire
};

map<string, int> ids;
vector<Reaction> reactions;
vector<long long> price;      // -1 daca nu se poate cumpara
int K;                        // numarul de molecule

// ---- starea simularii ----
vector<int> best;             // reactia aleasa pentru molecula (-1 = se cumpara)
vector<long long> inv;
vector<int> seq;
vector<char> inProg;          // protectie impotriva ciclurilor
long long totalCost, limitCost;
bool aborted;

int getId(const string& name){
    auto it = ids.find(name);
    if(it != ids.end()) return it->second;
    int id = ids.size();
    ids[name] = id;
    return id;
}

void add(Mols& list, int id, long long count){
    for(auto& p : list){
        if(p.first == id){
            p.second += count;
            return;
        }
    }
    list.push_back({id, count});
}

string readLine(){
    string str;
    while(getline(fin, str)){
        while(!str.empty() && (str.back() == '\r' || str.back() == ' ')) str.pop_back();
        if(!str.empty()) return str;
    }
    return "";
}

void parseReaction(const string& line, Reaction& r){
    r.text = line;
    string tmp = line;
    for(char& ch : tmp){
        if(ch == ',') ch = ' ';
    }
    stringstream rs(tmp);
    string word, name;
    bool afterArrow = false;
    while(rs >> word){
        if(word == "->"){
            afterArrow = true;
            continue;
        }
        rs >> name;
        if(afterArrow) add(r.out, getId(name), stoll(word));
        else add(r.in, getId(name), stoll(word));
    }
}

// ================= SIMULARE =================

void ensure(int m, long long q, int depth);

// executa o reactie o data; ce lipseste din inventar se cumpara automat
void run(int r, int depth){
    if(aborted) return;
    const Reaction& re = reactions[r];

    for(int pass = 0; pass < 3; pass++){
        for(auto& p : re.in) ensure(p.first, p.second, depth);
    }
    if(aborted) return;

    for(auto& p : re.in){
        if(inv[p.first] < p.second){
            if(price[p.first] < 0){ aborted = true; return; } // nu se poate cumpara
            totalCost += (p.second - inv[p.first]) * price[p.first];
            inv[p.first] = p.second;
        }
        inv[p.first] -= p.second;
    }
    for(auto& p : re.out) inv[p.first] += p.second;
    seq.push_back(r);
    if(totalCost >= limitCost || seq.size() >= MAXSEQ) aborted = true;
}

void ensure(int m, long long q, int depth){
    if(inProg[m]) return;
    inProg[m] = 1;
    while(!aborted && inv[m] < q){
        if(best[m] < 0 || depth > 300) break; // se va cumpara automat
        run(best[m], depth + 1);
    }
    inProg[m] = 0;
}

// intoarce costul, sau -1 daca planul e invalid / nu bate limita
long long simulate(int target, long long need, int bestT, long long lim){
    inv.assign(K, 0);
    inProg.assign(K, 0);
    seq.clear();
    totalCost = 0;
    limitCost = lim;
    aborted = false;
    while(!aborted && inv[target] < need) run(bestT, 0);
    if(aborted) return -1;
    return totalCost;
}

// ================= COSTURI UNITARE =================

// cost estimat pentru o bucata din "mol" produsa de reactia r
// alpha = cat din valoarea subproduselor se scade (0 = deloc)
double unitVal(const vector<double>& uc, int r, int idx, double alpha){
    const Reaction& re = reactions[r];
    double ic = 0;
    for(auto& p : re.in){
        if(uc[p.first] >= INF) return INF;
        ic += p.second * uc[p.first];
    }
    double credit = 0;
    if(alpha > 0){
        for(size_t j = 0; j < re.out.size(); j++){
            if((int)j == idx || re.outNet[j] <= 0) continue;
            if(uc[re.out[j].first] < INF) credit += re.outNet[j] * uc[re.out[j].first];
        }
    }
    double num = ic - alpha * credit;
    if(num < 0) num = 0;
    return num / re.outNet[idx];
}

int main(){
    auto globalStart = chrono::steady_clock::now();
    (void)globalStart;

    // 1. molecule cumparabile
    int m = stoi(readLine());
    vector<pair<int, long long>> buys;
    for(int i = 0; i < m; i++){
        stringstream ss(readLine());
        long long cost;
        string arrow, name;
        ss >> cost >> arrow >> name;
        buys.push_back({getId(name), cost});
    }

    // 2. reactii
    int r = stoi(readLine());
    reactions.resize(r);
    for(int i = 0; i < r; i++) parseReaction(readLine(), reactions[i]);

    // 3. tinte
    int n = stoi(readLine());
    vector<pair<int, long long>> targets;
    for(int i = 0; i < n; i++){
        stringstream ss(readLine());
        long long count;
        string name;
        ss >> count >> name;
        targets.push_back({getId(name), count});
    }

    K = ids.size();
    price.assign(K, -1);
    for(auto& b : buys) price[b.first] = b.second;

    // productia neta a fiecarei iesiri + lista de producatori pentru fiecare molecula
    // producers[mol] = (reactie, index in lista de iesiri)
    vector<vector<pair<int, int>>> producers(K);
    for(int i = 0; i < r; i++){
        Reaction& re = reactions[i];
        re.outNet.assign(re.out.size(), 0);
        for(size_t j = 0; j < re.out.size(); j++){
            long long consumed = 0;
            for(auto& p : re.in){
                if(p.first == re.out[j].first) consumed = p.second;
            }
            re.outNet[j] = re.out[j].second - consumed;
            if(re.outNet[j] > 0) producers[re.out[j].first].push_back({i, (int)j});
        }
    }

    // optiunile pentru fiecare molecula: -1 (cumpara) sau una din reactiile care o produc
    vector<vector<int>> options(K);
    for(int mol = 0; mol < K; mol++){
        if(price[mol] >= 0) options[mol].push_back(-1);
        for(auto& pr : producers[mol]) options[mol].push_back(pr.first);
    }

    // 4. mai multe puncte de start, cu valori diferite pentru subproduse
    const double alphas[3] = {0.0, 0.5, 1.0};
    vector<vector<double>> startCost(3);
    vector<vector<int>> startBest(3);
    for(int a = 0; a < 3; a++){
        vector<double>& uc = startCost[a];
        vector<int>& bs = startBest[a];
        uc.assign(K, INF);
        bs.assign(K, -1);
        for(int mol = 0; mol < K; mol++){
            if(price[mol] >= 0) uc[mol] = price[mol];
        }
        for(int iter = 0; iter < 300; iter++){
            bool changed = false;
            for(int mol = 0; mol < K; mol++){
                for(auto& pr : producers[mol]){
                    double val = unitVal(uc, pr.first, pr.second, alphas[a]);
                    if(val < uc[mol] - 1e-9){
                        uc[mol] = val;
                        bs[mol] = pr.first;
                        changed = true;
                    }
                }
            }
            if(!changed) break;
        }
    }

    // 5. fiecare tinta
    for(int t = 0; t < n; t++){
        int target = targets[t].first;
        long long need = targets[t].second;
        auto caseStart = chrono::steady_clock::now();
        auto elapsed = [&](){
            return chrono::duration<double>(chrono::steady_clock::now() - caseStart).count();
        };

        // cea mai buna pornire dintre cele 3
        long long cur = BIG;
        vector<int> bestStart;
        int bestStartT = -1;
        for(int a = 0; a < 3; a++){
            int bt = -1;
            double bv = INF;
            for(auto& pr : producers[target]){
                double v = unitVal(startCost[a], pr.first, pr.second, alphas[a]);
                if(v < bv){ bv = v; bt = pr.first; }
            }
            if(bt < 0) continue;
            best = startBest[a];
            long long c = simulate(target, need, bt, cur);
            if(c >= 0 && c < cur){
                cur = c;
                bestStart = startBest[a];
                bestStartT = bt;
            }
        }

        if(bestStartT < 0){
            // nu s-a gasit nimic valid; scriem un paragraf gol
            if(t + 1 < n) fout << "\n";
            cerr << "caz " << t + 1 << ": nu s-a gasit plan\n";
            continue;
        }

        best = bestStart;
        int bestT = bestStartT;

        // 6. cautare locala: schimba cate o alegere si pastreaza daca scade costul
        bool improved = true;
        while(improved && elapsed() < TIME_PER_CASE){
            improved = false;

            // alternativele pentru reactia tintei
            for(auto& pr : producers[target]){
                if(pr.first == bestT) continue;
                int old = bestT;
                bestT = pr.first;
                long long c = simulate(target, need, bestT, cur);
                if(c >= 0 && c < cur){ cur = c; improved = true; }
                else bestT = old;
            }

            // alternativele pentru celelalte molecule
            for(int mol = 0; mol < K && elapsed() < TIME_PER_CASE; mol++){
                if(options[mol].size() < 2) continue;
                for(int o : options[mol]){
                    if(o == best[mol]) continue;
                    int old = best[mol];
                    best[mol] = o;
                    long long c = simulate(target, need, bestT, cur);
                    if(c >= 0 && c < cur){ cur = c; improved = true; }
                    else best[mol] = old;
                }
            }
        }

        // regenereaza secventa finala
        long long fin_cost = simulate(target, need, bestT, BIG);
        for(int idx : seq) fout << reactions[idx].text << "\n";
        if(t + 1 < n) fout << "\n";
        cerr << "caz " << t + 1 << ": cost " << fin_cost << ", " << seq.size() << " reactii\n";
    }
}