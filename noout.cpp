#include <cassert>
#include <cstdio>
#include <ios>
#include <utility>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "splender.h"
#pragma GCC diagnostic pop
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

constexpr int deckn = 3, stage_n = 3, res_max = 3, reserved = 3, goal = 15;

// parameters
constexpr int depth[deckn]{15, 10, 7};
constexpr double rate[deckn]{0.96, 0.88, 0.7};
constexpr double buybase[stage_n] = {4, 5, 6}, resbase[stage_n] = {6, 7, 8},
                 part = 0.4,
                 dec_jk = 0.7; // for the second deck
constexpr double stage_parm[stage_n][3]{
{1.2, 0.4, 0.9}, {0.8, 1.1, 1}, {0.5, 1, 1}};
constexpr double res_parm[res_max]{0.75, 1, 1.25};
constexpr double majbuf = 3, dec_ratio = 0.75;

constexpr int TAKE = 1, BUY = 2, RES = 3;

namespace Gem {
constexpr int all = 6, normal = 5, num = 4, joker = 5;
constexpr int hold_max = 10;
int remain[Gem::all];
}; // namespace Gem

int stage, major = -1;
double demand[Gem::normal];
vector<pair<int, int>> resource[Gem::normal];

class Board {
  public:
    static constexpr int row = 4;
    card *operator[](int r) { return table[r]; }

  private:
    card table[deckn][row]; //,top[3][5];//0.9 0.7 0.5 0.3 0.1
};
Board table;

class Deck {
    vector<card> deck;
    int top_i;

  public:
    void init(vector<card> &orig) {
        deck = std::move(orig);
        top_i = 4;
    }
    card pop() {
        if (top_i >= static_cast<int>(deck.size())) {
            cerr << "deck empty\n";
            return (card){0, {0, 0, 0, 0, 0}, -1, 0};
        }
        return deck[top_i++];
    }
    card operator[](int i) { // 0 for top, not vector's 0
        i += top_i;
        if (i >= static_cast<int>(deck.size()) || i < top_i) {
            cerr << "invalid index\n";
            exit(1);
        }
        return deck[i];
    }
    bool empty() {
        return top_i == static_cast<int>(deck.size()) ? true : false;
    }
};
Deck deck[3];

struct ResCard {
    card c;
    int lv;
    ResCard(card c, int lv) : c(c), lv(lv) {}
};

struct profile {
    int gem[Gem::all], gem_sum;
    int bns[Gem::normal], bns_sum;
    int pts;
    int step[deckn][table.row];
    vector<ResCard> res;
    void init() {
        fill(gem, gem + Gem::all, 0);
        fill(bns, bns + Gem::normal, 0);
        bns_sum = 0;
        pts = 0;
        res.clear();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j) step[i][j] = 0;
    }
    void buy(card c) {
        for (int k = 0; k < Gem::normal; ++k) {
            int spent = max(0, c.cost[k] - bns[k]);
            Gem::remain[k] += min(spent, gem[k]);
            gem[k] -= spent;
            gem_sum -= spent;
            if (gem[k] < 0) {
                gem[Gem::joker] += gem[k];
                Gem::remain[Gem::joker] -= gem[k];
                gem[k] = 0;
            }
        }
        pts += c.score;
        ++bns[c.gem];
        ++bns_sum;
    }
    void resv(card c, int lv) {
        res.push_back(ResCard(c, lv));
        ++gem[Gem::joker];
        ++gem_sum;
        --Gem::remain[Gem::joker];
    }
};
profile my, op;

double fac[deckn][depth[0]];
void fac_init() {
    assert(depth[0] >= depth[1] && depth[0] >= depth[2]);
    for (int lv = 0; lv < deckn; ++lv) {
        for (int i = 0; i <= table.row; ++i) fac[lv][i] = 1; // first five
        for (int i = table.row + 1; i < depth[lv]; ++i)
            fac[lv][i] = fac[lv][i - 1] * rate[lv];
    }
}

bool reservable(profile &p) {
    return p.res.size() < res_max && p.gem_sum < 10 &&
           Gem::remain[Gem::joker] > 0;
}

// cal_step
struct Step {
    bool toofar;
    int pickn;
    int gem[Gem::normal];
    Step() : toofar(false), pickn(0), gem() {}
};

Step cal_step(const card c,
              const profile &p) { // problem: no joker and no bound
    Step ret;
    // check whether toofar
    for (int i = 0; i < Gem::normal; ++i) {
        /*
        if(c.cost[i]>=4){
                int d=c.cost[i]-p.gem[i]-p.bns[i]-p.gem[Gem::joker];
                if(d<=3)
                if(d<=3&&Gem::remain>=d&&reservable(p))

        }
        */
        if ((c.cost[i] == 5 && p.bns[i] < 1) ||
            (c.cost[i] == 6 && p.bns[i] < 2) ||
            (c.cost[i] == 7 && p.bns[i] < 3)) {
            ret.toofar = true;
            return ret;
        }
    }

    // pick
    int diff[Gem::normal]{}, match = 0, color_match = 0;
    if (my.gem_sum <= 7) ret.pickn = 3;
    else if (my.gem_sum == 8) ret.pickn = 2;

    for (int i = 0; i < Gem::normal; ++i) {
        diff[i] = c.cost[i] - p.gem[i] - p.bns[i];
        if (diff[i] > 0) {
            int nmatch = min(diff[i], Gem::remain[i]);
            if (ret.pickn == 2 && match < 2 && nmatch > match) {
                ret.gem[color_match] = 0;
                ret.gem[i] = min(2, nmatch);
                match = ret.gem[i];
                color_match = i;
            }
            else if (ret.pickn == 3 && match < 3 && nmatch > 0) {
                ++match;
                color_match = i;
                ++ret.gem[i];
            }
        }
        else if (diff[i] < 0) diff[i] = 0;
    }

    if (ret.pickn == 3 && match == 1) {
        if (diff[color_match] >= 2 && Gem::remain[color_match] >= 2) {
            ret.pickn = 2;
            ret.gem[color_match] = 2;
        }
    }
    else if (match == 0 || (ret.pickn == 2 && ret.gem[color_match] < 2)) {
        ret.pickn = 0;
        /*
for (int i = 0; i < Gem::normal; ++i)
    ret.gem[i] = min(diff[i], 1);
        */
    }

    return ret;
}

bool iskeycard(card c) {
    int m = 0;
    for (int i = 1; i < Gem::normal; ++i) {
        if (c.cost[i] > c.cost[m]) m = i;
    }
    if (m == major && c.cost[m] >= 6) return true;
    return false;
}

double score1(const card &cd) {
    static constexpr double c[6] = {0,       1 * 1,    1.2 * 2,
                                    1.3 * 3, 1.45 * 4, 1.55 * 5};
    double deduction = 0;
    // cout << cd.gem << ": ";
    for (int i = 0; i < Gem::normal; ++i) {
        // cout << cd.cost[i] << ' ';
        int cost = cd.cost[i] - my.bns[i];
        if (cost > 0) deduction += c[cost];
    }
    // cout << 10 - score << endl;
    return 10 - deduction;
}

double score2(const card &cd) {
    static constexpr double c1 = 1.2, c2 = 1.3; // c1*n - c2*pt
    double deduction = 0;
    // cout << cd.gem << ": ";
    for (int i = 0; i < Gem::normal; ++i) {
        int cost = cd.cost[i] - my.bns[i];
        // cout << cost<< ' ';
        if (cost > 0) deduction += c1 * cost;
    }
    // cout << score - c2 * cd.score << '\n';
    return 10 - (deduction - c2 * cd.score);
}

double score3(const card &cd) {
    static constexpr double c1 = 1.2, c2 = 1.3;
    double deduction = 0;
    // cout << cd.gem << ": ";
    for (int i = 0; i < Gem::normal; ++i) {
        int cost = cd.cost[i] - my.bns[i];
        // cout << cost << ' ';
        if (cost > 0) deduction += c1 * (cost - 1);
    }
    // cout << deduction - c2 * cd.score << '\n';
    return 10 - (deduction - c2 * cd.score);
}

double (*scoref[])(const card &) = {score1, score2, score3};

void print_top(vector<card> &st, int lv) {
    for (int t = 0; t < Gem::normal; ++t) {
        cout << t << ": ";
        for (int i = 0; i < depth[lv]; ++i) {
            if (st[i].gem == t) cout << i << ' ';
        }
        cout << endl;
    }
}

void init(vector<card> stack_1, vector<card> stack_2, vector<card> stack_3) {
    // cout << setprecision(2) << fixed;
    //  init var
    fac_init();
    my.init();
    op.init();
    stage = 0;
    fill(demand, demand + 5, 0);

    // analyze here because index of deck start from 4
    // supply 1
    double supply[Gem::normal]{};
    for (int i = 0; i < depth[0]; ++i) {
        supply[stack_1[i].gem] += fac[0][i] * score1(stack_1[i]);
        // cout << stack_1[i].gem << " += " << score1(stack_1[i])
        //<< '\n';
    }

    // cout << "Lv 0" << endl;
    // print_top(stack_1, 0);

    double sup1[Gem::normal] = {supply[0], supply[1], supply[2], supply[3],
                                supply[4]};
    // supply and demand 2
    for (int i = 0; i < depth[1]; ++i) {
        double v = fac[1][i] * part * score2(stack_2[i]);
        supply[stack_2[i].gem] += v;
        for (int type = 0; type < Gem::normal; ++type) {
            if (stack_2[i].cost[type] > 3)
                demand[type] += v * stack_2[i].cost[type];
        }
        // cout << stack_2[i].gem << " += " << score2(stack_2[i])
        //<< '\n';
    }

    /*
cout << "Supply:\n";
for (int i = 0; i < Gem::normal; ++i) {
    cout << i << "  Lv 0: " << sup1[i] << "\tLv 1: " << supply[i] - sup1[i]
         << "\tTotal: " << supply[i] << '\n';
}

cout << "demand before:\n";
for (double i : demand) cout << i << ' ';
cout << endl;
    */
    // demand 3
    for (int i = 0; i < depth[2]; ++i) {
        double v = fac[2][i] * score3(stack_3[i]);
        for (int type = 0; type < Gem::normal; ++type) {
            /*
cout << v *
        (stack_3[i].cost[type] < 3 ? 0 : stack_3[i].cost[type] - 3)
     << ' ';
                     */
            if (stack_3[i].cost[type] > 3)
                demand[type] += v * (stack_3[i].cost[type] - 3);
        }
        // cout << endl << endl;
    }
    // cout << "demand after:\n";
    // for (double i : demand) cout << i << ' ';
    // cout << endl;

    // for (int i = 0; i < Gem::normal; ++i) demand[i] *= supply[i];
    // cout << "demand times supply:\n";
    // for (double i : demand) cout << i << ' ';
    // cout << endl;
    //  decide major and minor by supply and demand
    major = 0;
    for (int i = 1; i < Gem::normal; ++i) {
        if (demand[i] > demand[major]) {
            major = i;
        }
    }
    // cout << "Major: " << major << endl;
    for (int i = 0; i < Gem::normal; ++i) {
        if (i == major) demand[i] = majbuf;
        else demand[i] = 1;
    }

    // set table and deck
    fill(Gem::remain, Gem::remain + Gem::all, 4);
    for (int i = 0; i < table.row; ++i) {
        table[0][i] = stack_1[i];
        table[1][i] = stack_2[i];
        table[2][i] = stack_3[i];
    }
    deck[0].init(stack_1); // stack "moved"
    deck[1].init(stack_2);
    deck[2].init(stack_3);
}

pair<int, int> get_index(profile &p, int cid) {
    for (int i = 0; i < deckn; ++i) {
        for (int j = 0; j < table.row; ++j) {
            if (table[i][j].id == cid) {
                return pair<int, int>(i, j);
            }
        }
    }
    for (int i = 0; i < p.res.size(); ++i) {
        if (p.res[i].c.id == cid) return pair<int, int>(reserved, i);
    }
    cerr << "cannot get index\n";
    exit(1);
}
bool iscard(card c) { return c.gem != -1; }

void update(profile &p, pair<int, int> card_index, const int type) {
    if (card_index.first == reserved) { // buy reserved card
        assert(type == BUY);
        auto it = p.res.begin() + card_index.second;
        p.buy(it->c);
        p.res.erase(it);
    }
    else {
        const int i = card_index.first, j = card_index.second;
        if (type == BUY) {
            p.buy(table[i][j]);
            table[i][j] = deck[i].pop();
        }
        else { // type==RES
            p.resv(table[i][j], i);
            table[i][j] = deck[i].pop();
        }
    }
}

int lack(profile &p, card c, int joker_max) {
    // cout << "[";
    // for (int i = 0; i < Gem::normal; ++i) cout << c.cost[i] << ' ';
    // cout << "]\n";
    int jk = 0;
    for (int i = 0; i < Gem::normal; ++i) {
        if (c.cost[i] > p.gem[i] + p.bns[i]) {
            jk += c.cost[i] - p.gem[i] - p.bns[i];
        }
    }
    if (jk > joker_max) {
        // cout << "cannot afford\n";
        return jk - joker_max;
    }
    return 0;
}

struct move buy(pair<int, int> card_index) {
    card c;
    if (card_index.first < reserved) {
        c = table[card_index.first][card_index.second];
    }
    else {
        c = my.res[card_index.second].c;
    }
    update(my, card_index, BUY);
    demand[c.gem] *= dec_ratio;
    return (struct move){BUY, c.id, {0, 0, 0, 0, 0}};
}

struct move reserve(pair<int, int> card_index) {
    card c = table[card_index.first][card_index.second];
    update(my, card_index, RES);
    return (struct move){RES, c.id, {0, 0, 0, 0, 0}};
}

struct move take(int g[]) {
    for (int i = 0; i < Gem::normal; ++i) {
        if (g[i] != 0) {
            my.gem[i] += g[i];
            my.gem_sum += g[i];
            Gem::remain[i] -= g[i];
        }
    }
    return {TAKE, 0, {g[0], g[1], g[2], g[3], g[4]}};
}

void print(profile &p) {
    cout << "gem: ";
    for (int i = 0; i < Gem::all; ++i) cout << p.gem[i] << ' ';
    cout << endl;
    cout << "bonus: ";
    for (int i = 0; i < Gem::normal; ++i) cout << p.bns[i] << ' ';
    cout << endl;
    cout << "points: " << p.pts << endl;
    cout << "card reserved:\n";
    for (auto cd : p.res) {
        cout << "Gem: " << cd.c.gem << " Cost: ";
        for (auto i : cd.c.cost) cout << i << ' ';
        cout << endl;
    }
    cout << endl;
}

void print_table() {
    for (int i = 0; i < deckn; ++i) {
        for (int j = 0; j < table.row; ++j) {
            cout << table[i][j].id << " [";
            for (int k = 0; k < Gem::normal; ++k) {
                cout << table[i][j].cost[k];
                if (k != Gem::normal - 1) cout << ' ';
            }
            cout << "][";
            cout << table[i][j].gem << ' ' << table[i][j].score << "] ";
        }
        cout << endl;
    }
    cout << endl;
}

bool allzero(int a[]) {
    for (int i = 0; i < Gem::normal; ++i)
        if (a[i] != 0) return false;
    return true;
}

void print_step(Step step) {
    cout << "step.toofar = " << ((step.toofar) ? "true" : "false")
         << " step.pickn = " << step.pickn << endl;
    for (int i = 0; i < Gem::normal; ++i) cout << step.gem[i] << ' ';
    cout << endl;
}

struct move player_move(struct move m) {
    // record opponent's move
    if (m.type != 0) { // eupdate board
        switch (m.type) {
        case TAKE:
            for (int i = 0; i < Gem::normal; ++i) {
                op.gem[i] += m.gem[i];
                Gem::remain[i] -= m.gem[i];
            }
            break;
        case BUY:
            update(op, get_index(op, m.card_id), BUY);
            break;
        case RES:
            update(op, get_index(op, m.card_id), RES);
            break;
        }
    }
    // cout << "Opponent's:\n";
    // print(op);

    // cout << "Mine:\n";
    // print(my);
    // print_table();

    // compute value and step
    vector<pair<int, int>> card_index;

    // change stage
    if (stage == 0 && my.bns[major] >= 3) stage = 1;
    else if (my.pts > 8) stage = 2;

    // reserved cards
    double value[deckn + 1][table.row]{};

    for (int i = 0; i < my.res.size(); ++i) {
        value[reserved][i] = scoref[my.res[i].lv](my.res[i].c) *
                             stage_parm[stage][my.res[i].lv] *
                             res_parm[my.res.size() - 1];
        if (my.res[i].lv == 2 && iskeycard(my.res[i].c))
            value[reserved][i] *= 3;
        else if (my.res[i].lv <= 1)
            value[reserved][i] *= demand[my.res[i].c.gem];
        /*
if (my.res[i].lv == 0)
    value[reserved][i] = score1(my.res[i].c);
else if (my.res[i].lv == 1)
    value[reserved][i] = score2(my.res[i].c);
else value[reserved][i] = score3(my.res[i].c);
        */
        card_index.push_back(pair<int, int>(reserved, i));
    }

    for (int i = 0; i < Gem::normal; ++i) resource[i].clear();

    // table cards
    for (int lv = 0; lv < 2; ++lv) {
        for (int i = 0; i < table.row; ++i) {
            if (iscard(table[lv][i])) {
                value[lv][i] = scoref[lv](table[lv][i]) *
                               stage_parm[stage][lv] * demand[table[lv][i].gem];
                card_index.push_back(pair<int, int>(lv, i));
                resource[table[lv][i].gem].push_back(pair<int, int>(lv, i));
            }
        }
    }
    // lv = 2
    for (int i = 0; i < table.row; ++i) {
        if (iscard(table[2][i])) {
            value[2][i] = scoref[2](table[2][i]) * stage_parm[stage][2] *
                          (iskeycard(table[2][i]) ? 3 : 1);
            card_index.push_back(pair<int, int>(2, i));
            resource[table[2][i].gem].push_back(pair<int, int>(2, i));
        }
    }
    auto cmp = [value](pair<int, int> a, pair<int, int> b) -> bool {
        return value[a.first][a.second] > value[b.first][b.second];
    };
    for (int i = 0; i < Gem::normal; ++i) {
        sort(resource[i].begin(), resource[i].end(),
             cmp); // capture not sure
    }
    sort(card_index.begin(), card_index.end(),
         cmp); // capture not sure

    // choose best move

    // pair<int, int> bestres;
    int totake[5]{};
    int totake_n = 0, filled_n = 0;
    double dependency[5]{};
    int jkcnt[deckn + 1][table.row]{};
    // auto it = card_index.begin(); //fail after insert
    int cur = 0;
    for (double dep_set_n = 0, val;
         cur != card_index.size() &&
         value[card_index[cur].first][card_index[cur].second] > buybase[stage];
         ++cur) {
        auto it = card_index.begin() + cur;
        val = value[it->first][it->second];
        card c = (it->first == reserved) ? my.res[it->second].c
                                         : table[it->first][it->second];
        int dgem = lack(my, c, jkcnt[it->first][it->second]);
        if (dgem == 0) return buy(*it);
        else {
            if (it->first != reserved && reservable(my) &&
                val > resbase[stage]) {
                return reserve(*it);
            }
            else {
                Step step = cal_step(c, my);
                print_step(step);

                if (!step.toofar) {
                    if (totake_n == 0 && step.pickn != 0) {
                        totake_n = step.pickn;
                        for (int i = 0; i < Gem::normal; ++i) {
                            if (step.gem[i] != 0) {
                                totake[i] = step.gem[i];
                                filled_n += step.gem[i];
                            }
                        }
                    }
                    else if (!allzero(step.gem)) {
                        for (int i = 0; i < Gem::normal; ++i) {
                            dependency[i] += (double)step.gem[i] / dep_set_n;
                        }
                        ++dep_set_n;
                    }

                    if (jkcnt[it->first][it->second] < my.gem[Gem::joker]) {
                        // reshuffle by value
                        ++jkcnt[it->first][it->second];
                        value[it->first][it->first] -= dec_jk;
                        auto pos = lower_bound(card_index.begin(),
                                               card_index.end(), *it, cmp);
                        card_index.insert(pos, *it);
                    }
                }
            }
        }
    }

    if (totake_n == 0 || filled_n < totake_n) { // only when
                                                // totake_n=3
        if (my.gem_sum == 8) {
            int a[Gem::normal] = {0, 1, 2, 3, 4};
            sort(a, a + Gem::normal, [dependency](int a, int b) {
                if (Gem::remain[a] < 2) return false;
                if (Gem::remain[b] < 2) return true;
                return dependency[a] > dependency[b];
            });
            if (Gem::remain[a[0]] >= 2) { // if first not ok, no one
                                          // ok
                totake[a[0]] = 2;
                return take(totake);
            }
        }
        else if (my.gem_sum < 8) {
            int a[Gem::normal] = {0, 1, 2, 3, 4};
            sort(a, a + Gem::normal, [dependency](int a, int b) {
                if (Gem::remain[a] == 0) return false;
                if (Gem::remain[b] == 0) return true;
                return dependency[a] > dependency[b];
            });
            for (int i = 0;
                 i < Gem::normal && Gem::remain[i] > 0 && filled_n < 3; ++i) {
                if (totake[a[i]] == 0) {
                    totake[a[i]] = 1;
                    ++filled_n;
                }
            }
            if (filled_n == 3) return take(totake);
            else {
                fill(totake, totake + Gem::normal, 0);
                for (int i = 0; i < 2; ++i) {
                    if (Gem::remain[a[i]] >= 2) {
                        totake[a[i]] = 2;
                        return take(totake);
                    }
                }
            }
        }
    }
    else {
        return take(totake);
    }
    for (auto iter = card_index.begin() + cur; iter != card_index.end();
         ++iter) {
        card c;
        if (iter->first != reserved) {
            c = table[iter->first][iter->second];
        }
        else {
            c = my.res[iter->second].c;
        }
        if (lack(my, c, jkcnt[iter->first][iter->second]) == 0)
            return buy(*iter);
    }
    // cannot buy anything try reserve
    if (reservable(my)) {
        for (auto iter = card_index.begin() + cur; iter != card_index.end();
             ++iter) {
            if (iter->first != reserved) {
                return reserve(*iter);
            }
        }
    }

    // no valid operation
    cerr << "No valid operation\n";
    m.type = 3;
    m.card_id = 89;
    return m;
}