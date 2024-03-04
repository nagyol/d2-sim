#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <functional>
#include <cassert>

using namespace std;

//Namespace describition auxiliary utilities - pretty-printers, RNGs, structure explorers etc. etc.
namespace utils
{
    mt19937_64 get_RNG()
    {
        random_device rd;
        mt19937_64 g(rd());
        return g;
    };

    //Get a pair of two distinct numbers chosen uar from [1, upper_bound].
    pair<int, int> get_random_pair(const int upper_bound, mt19937_64 &RNG)
    {
        uniform_int_distribution<> distrib1(1, upper_bound);
        uniform_int_distribution<> distrib2(1, upper_bound - 1);
        int first = distrib1(RNG);
        int second = distrib2(RNG);
        if (second >= first)
            ++second;
        return {first, second};
    }

    //Get a number chosen uar from [1, upper_bound].
    int get_random_number(int upper_bound, mt19937_64 &RNG)
    {
        uniform_int_distribution<> distrib(1, upper_bound);
        return distrib(RNG);
    }

    //Returns a vector containing elements of permutation that are on the same cycle as elem.
    //Cycles start with elem.
    vector<int> get_cycle_of_element(const vector<int> &perm, const int elem)
    {
        vector<int> cycle;
        int first = elem;
        int current = perm[elem - 1];
        cycle.push_back(first);
        while (first != current)
        {
            cycle.push_back(current);
            current = perm[current - 1];
        }
        return (cycle);
    }

    //Returns of a vector of all cycles. Cycles themselves in the same form as in get_cycle_of_element(...).
    vector<vector<int>> get_cycles(const vector<int> &perm)
    {
        vector<vector<int>> cycles;

        vector<int> indices(perm.size());
        iota(begin(indices), end(indices), 0);

        for (size_t i = 0; i < perm.size(); ++i)
        {
            if (indices[i] == -1)
                continue;
            vector<int> cycle;
            int first = i + 1;
            int current = perm[i];
            cycle.push_back(first);
            indices[first - 1] = -1;
            while (first != current)
            {
                cycle.push_back(current);
                indices[current - 1] = -1;
                current = perm[current - 1];
            }
            cycles.push_back(cycle);
        };
        return cycles;
    };

    vector<vector<int>> get_cycles_sorted(const vector<int> &perm)
    {
        auto cycles = get_cycles(perm);
        sort(cycles.begin(), cycles.end(), [](const vector<int> & a, const vector<int> & b){ return a.size() > b.size(); });
        return cycles;
    }

    //[Helper] Print permutation element by element.
    void print_permutation(const vector<int> &perm)
    {
        for (auto elem : perm)
            cout << elem << ", ";
        cout << "\n";
    }

    //[Helper] Print distribution element by element.
    void print_distribution(const vector<double> &dist)
    {
        for (auto elem : dist)
            cout << elem << ", ";
        cout << "\n";
    }

    //[Helper] Pretty-printer for cycles and elements withing them.
    //Example output: (1,2,3,)(4,) . If parsing, watch out for trailing commas.
    void print_cycles(const vector<vector<int>> &cycles)
    {
        for (auto &cycle : cycles)
        {
            cout << "(";
            for (auto elem : cycle)
            {
                cout << elem << ",";
            }
            cout << ")";
        }
        cout << "\n";
    }

    //Generate a vector of all transposition constructable on {1, ..., graph_size}.
    vector<pair<int, int>> get_all_transpositions(int graph_size)
    {
        vector<pair<int, int>> transpositions;
        for (int i = 1; i <= graph_size; ++i)
        {
            for (int j = graph_size; j > i; --j)
            {
                pair<int, int> transp = {i, j};
                transpositions.push_back(transp);
            }
        }
        return transpositions;
    }

    //Return TV distance from Unif(1, ..., dist.size()), calculated as a half of L_1 distance.
    double get_TV_distance_from_uniform(vector<double> dist)
    {
        double TVD = 0.;
        for (auto i : dist)
        {
            TVD += abs(i - 1. / dist.size());
        }
        return TVD / 2;
    }

}

////////////////////////////////////////////////////

//Namespace describing different generators of the initial graph.
namespace init_graph
{
    //Generate the initial graph as a random permutation.
    vector<int> generate_random_permutation(const int size, mt19937_64 &RNG)
    {
        vector<int> perm;
        perm.reserve(size);
        for (int i = 1; i <= size; ++i)
            perm.push_back(i);
        shuffle(perm.begin(), perm.end(), RNG);
        return perm;
    };

    //Generate an initial graph with size vertices that consists solely of singletons.
    vector<int> generate_singletons(const int size)
    {
        vector<int> sing;
        sing.reserve(size);
        for (int i = 1; i <= size; ++i)
            sing.push_back(i);
        return sing;
    };
}

////////////////////////////////////////////////////

//Namespace describing functions providing the evolution of both the graph and RW.
namespace dynamics
{
    //Apply permutation (first, second) on &perm.
    void apply_transposition(const int first, const int second, vector<int> &perm)
    {
        int temp = perm[first - 1];
        perm[first - 1] = perm[second - 1];
        perm[second - 1] = temp;
    };

    //Initialize the Infinite Speed Random Walk. Distribution is uniform over the cycle that contains init_pos.
    vector<double> init_ISRW(const vector<int> &perm, const int init_pos)
    {
        vector<double> distribution;
        distribution.reserve(perm.size());
        for (size_t i = 0; i < perm.size(); ++i)
            distribution.push_back(0);
        auto init_cycle = utils::get_cycle_of_element(perm, init_pos);
        for (auto c_elems : init_cycle)
            distribution[c_elems - 1] = 1. / init_cycle.size();
        return distribution;
    };

    //Evolves the distribution of the ISRW if the dynamics of the graph can be described by a transposition.
    //Notice that for transposition-generated graph dynamics we can easily see which entries in the distribution will change.
    void evolve_ISRW_transp(vector<int> &perm, vector<double> &dist, const int trans1, const int trans2)
    {
        auto pre_cycle1 = utils::get_cycle_of_element(perm, trans1);
        auto pre_cycle2 = utils::get_cycle_of_element(perm, trans2);
        apply_transposition(trans1, trans2, perm);
        if (find(pre_cycle1.begin(), pre_cycle1.end(), trans2) == pre_cycle1.end())
        {
            double sum = dist[trans1 - 1] * pre_cycle1.size() + dist[trans2 - 1] * pre_cycle2.size();
            auto post_cycle = utils::get_cycle_of_element(perm, trans1);
            for (auto elem : post_cycle)
                dist[elem - 1] = sum / post_cycle.size();
        }
    };

}

////////////////////////////////////////////////////

//Namespace concerned with functions that carry out the concrete simulation scenarios.
//Cf. namespace dynamics that contains general transformations of RWs or graphs.
namespace simulation
{
    struct properties
    {
        size_t steps;
        double epsilon;
        string sep = " ";
        vector<int> init_graph;
        bool DEBUG = false;
        bool print_cycles = false;
        mt19937_64 RNG;
        bool fragmentation = true;
        bool coag_UOC = false;
    };

    pair<int, int> get_UOCcoagulative_transposition(const vector<int> &graph, simulation::properties &conf)
    {
        auto cycles = utils::get_cycles(graph);
        auto to_merge = utils::get_random_pair(cycles.size(), conf.RNG);
        return pair<int, int>(cycles.at(to_merge.first - 1).at(0), cycles.at(to_merge.second - 1).at(0));
    }

    pair<int, int> get_UOVcoagulative_transposition(const vector<int> &graph, simulation::properties &conf)
    {
        auto all_cycles = utils::get_cycles(graph);
        auto to_merge = utils::get_random_pair(graph.size(), conf.RNG);
        auto first_cycle = utils::get_cycle_of_element(graph, to_merge.first);
        if (find(first_cycle.begin(), first_cycle.end(), to_merge.second) != first_cycle.end())
        {
            if (conf.DEBUG)
                cout << "#CONFLICT!\n";

            vector<int> flattened_cycles;
            flattened_cycles.reserve(graph.size());
            for (size_t cycle_num = 0; cycle_num < all_cycles.size(); ++cycle_num)
            {
                flattened_cycles.insert(flattened_cycles.end(), all_cycles.at(cycle_num).size(), cycle_num);
                if (conf.DEBUG)
                    cout << "size Cycle" << cycle_num << ": " << all_cycles.at(cycle_num).size() << "\n";
            }
            if (conf.DEBUG)
                cout << flattened_cycles.size() << "  " << all_cycles.size() << "\n";
            auto first = utils::get_random_number(graph.size(), conf.RNG) - 1;
            auto second = utils::get_random_number(graph.size() - all_cycles.at(flattened_cycles.at(first)).size(), conf.RNG) - 1;
            if (flattened_cycles.at(second) >= flattened_cycles.at(first))
            {
                if (conf.DEBUG)
                    cout << "#CONFLICT 2! " << first << "  " << second << "\n";
                second = second + all_cycles.at(flattened_cycles.at(first)).size();
            }
            return {all_cycles.at(flattened_cycles.at(first)).at(0), all_cycles.at(flattened_cycles.at(second)).at(0)};
        }

        return to_merge;
    }

    pair<int, int> get_GD_generator(properties &conf, vector<int> &graph)
    {
        if (conf.fragmentation)
        {
            return utils::get_random_pair(graph.size(), conf.RNG);
        }
        else if (conf.coag_UOC)
        {
            return get_UOCcoagulative_transposition(graph, conf);
        }
        else
        {
            return get_UOVcoagulative_transposition(graph, conf);
        };
    }

    void simulate_quenched(simulation::properties &conf)
    {
        auto graph_size = conf.init_graph.size();
        auto graph = conf.init_graph;
        auto dist = dynamics::init_ISRW(graph, utils::get_random_pair(graph_size, conf.RNG).first);
        double TVD = utils::get_TV_distance_from_uniform(dist);

        if (conf.DEBUG)
        {
            utils::print_cycles(utils::get_cycles(graph));
            utils::print_distribution(dist);
        };
        cout << "0" << conf.sep
             << TVD << conf.sep
             << (conf.print_cycles ? to_string(utils::get_cycles(graph).size()) : "") << "\n"
             << (conf.DEBUG ? "\n" : "");

        auto transp = get_GD_generator(conf, graph);

        for (size_t cur_step = 1; cur_step <= conf.steps; ++cur_step)
        {
            if (TVD <= conf.epsilon)
                break;

            transp = get_GD_generator(conf, graph);

            if (conf.DEBUG)
                cout << "TSP: " << transp.first << " " << transp.second << "\n";

            dynamics::evolve_ISRW_transp(graph, dist, transp.first, transp.second);
            TVD = utils::get_TV_distance_from_uniform(dist);

            if (conf.DEBUG)
            {
                utils::print_cycles(utils::get_cycles(graph));
                utils::print_distribution(dist);
            };
            cout << cur_step << conf.sep
                 << TVD << conf.sep
                 << (conf.print_cycles ? to_string(utils::get_cycles(graph).size()) : "") << "\n"
                 << (conf.DEBUG ? "\n" : "");
        }
    }

    void simulate_quenched_biggest_ATM(simulation::properties &conf)
    {
        auto graph_size = conf.init_graph.size();
        auto graph = conf.init_graph;
        auto dist = dynamics::init_ISRW(graph, utils::get_random_pair(graph_size, conf.RNG).first);
        double TVD = utils::get_TV_distance_from_uniform(dist);

        if (conf.DEBUG)
        {
            utils::print_cycles(utils::get_cycles(graph));
            utils::print_distribution(dist);
        };
        cout << "0" << conf.sep
             << TVD << conf.sep
             << (conf.print_cycles ? to_string(utils::get_cycles(graph).size()) : "") << "\n"
             << (conf.DEBUG ? "\n" : "");

        auto transp = get_GD_generator(conf, graph);

        for (size_t cur_step = 1; cur_step <= conf.steps; ++cur_step)
        {
            if (TVD <= conf.epsilon)
                break;

            transp = get_GD_generator(conf, graph);

            if (conf.DEBUG)
                cout << "TSP: " << transp.first << " " << transp.second << "\n";

            dynamics::evolve_ISRW_transp(graph, dist, transp.first, transp.second);

            dist = dynamics::init_ISRW(graph, utils::get_cycles_sorted(graph)[0][0]);
            TVD = utils::get_TV_distance_from_uniform(dist);

            if (conf.DEBUG)
            {
                utils::print_cycles(utils::get_cycles(graph));
                utils::print_distribution(dist);
            };
            cout << cur_step << conf.sep
                 << TVD << conf.sep
                 << (conf.print_cycles ? to_string(utils::get_cycles(graph).size()) : "") << "\n"
                 << (conf.DEBUG ? "\n" : "");
        }
    }

    void simulate_quenched_biggest_sofar(simulation::properties &conf)
    {
        auto graph_size = conf.init_graph.size();
        auto graph = conf.init_graph;
        auto dist = dynamics::init_ISRW(graph, utils::get_random_pair(graph_size, conf.RNG).first);
        double TVD = utils::get_TV_distance_from_uniform(dist);

        if (conf.DEBUG)
        {
            utils::print_cycles(utils::get_cycles(graph));
            utils::print_distribution(dist);
        };
        cout << "0" << conf.sep
             << TVD << conf.sep
             << (conf.print_cycles ? to_string(utils::get_cycles(graph).size()) : "") << "\n"
             << (conf.DEBUG ? "\n" : "");

        auto transp = get_GD_generator(conf, graph);

        for (size_t cur_step = 1; cur_step <= conf.steps; ++cur_step)
        {
            if (TVD <= conf.epsilon)
                break;

            transp = get_GD_generator(conf, graph);

            if (conf.DEBUG)
                cout << "TSP: " << transp.first << " " << transp.second << "\n";

            dynamics::evolve_ISRW_transp(graph, dist, transp.first, transp.second);

            dist = dynamics::init_ISRW(graph, utils::get_cycles_sorted(graph)[0][0]);
            auto newTVD = utils::get_TV_distance_from_uniform(dist);
            if (newTVD < TVD) 
                TVD = newTVD;

            if (conf.DEBUG)
            {
                utils::print_cycles(utils::get_cycles(graph));
                utils::print_distribution(dist);
            };
            cout << cur_step << conf.sep
                 << TVD << conf.sep
                 << (conf.print_cycles ? to_string(utils::get_cycles(graph).size()) : "") << "\n"
                 << (conf.DEBUG ? "\n" : "");
        }
    }

    vector<double> simulate_quenched_r(simulation::properties &conf)
    {
        auto graph_size = conf.init_graph.size();
        auto graph = conf.init_graph;
        auto dist = dynamics::init_ISRW(graph, utils::get_random_pair(graph_size, conf.RNG).first);
        double TVD = utils::get_TV_distance_from_uniform(dist);
        vector<double> TVD_evolution;
        TVD_evolution.reserve(conf.steps + 1);

        if (conf.DEBUG)
        {
            utils::print_cycles(utils::get_cycles(graph));
            utils::print_distribution(dist);
        };

        TVD_evolution.push_back(TVD);

        auto transp = get_GD_generator(conf, graph);

        for (size_t cur_step = 1; cur_step <= conf.steps; ++cur_step)
        {
            if (TVD <= conf.epsilon)
                break;

            transp = get_GD_generator(conf, graph);

            if (conf.DEBUG)
                cout << "TSP: " << transp.first << " " << transp.second << "\n";

            dynamics::evolve_ISRW_transp(graph, dist, transp.first, transp.second);
            TVD = utils::get_TV_distance_from_uniform(dist);

            if (conf.DEBUG)
            {
                utils::print_cycles(utils::get_cycles(graph));
                utils::print_distribution(dist);
            };

            TVD_evolution.push_back(TVD);
        }
        return TVD_evolution;
    }

    void get_averaged_TVD_quenched(simulation::properties &CS, int runs)
    {
        vector<double> final_TVD(CS.steps + 1, 0.); //initial pos. + n steps
        for (int i = 1; i <= runs; ++i)
        {
            vector<double> TVD = simulate_quenched_r(CS);
            for (size_t j = 0; j < TVD.size(); ++j)
            {
                final_TVD.at(j) = ((final_TVD.at(j) * (i - 1)) + TVD.at(j)) / i;
            }
        cerr << "Finished run " << i << "\n"; 
        }
        for (size_t j = 0; j < CS.steps; ++j)
        {
            cout << j << CS.sep << final_TVD.at(j) << "\n";
            if (final_TVD.at(j) == 0.)
                break;
        }
    }

}

////////////////////////////////////////////////////

//Methods that test the program in simple and easy to follow situations.
namespace test
{
    void test_evolution()
    {
        vector<int> test = {4, 2, 5, 1, 3};
        //cycles: (1,4) (2) (3,5)
        auto dist = dynamics::init_ISRW(test, 1);
        //dist: 0.5, 0, 0, 0.5, 0
        utils::print_distribution(dist);
        cout << "1:" << utils::get_TV_distance_from_uniform(dist) << "\n";
        dynamics::evolve_ISRW_transp(test, dist, 4, 2);
        cout << "2:" << utils::get_TV_distance_from_uniform(dist) << "\n";
        //dist 1/3, 1/3, 0, 1/3, 0
        utils::print_distribution(dist);
        dynamics::evolve_ISRW_transp(test, dist, 4, 2);
        cout << "3:" << utils::get_TV_distance_from_uniform(dist) << "\n";
        //dist 1/3, 1/3, 0, 1/3, 0
        utils::print_distribution(dist);
        dynamics::evolve_ISRW_transp(test, dist, 2, 3);
        cout << "4:" << utils::get_TV_distance_from_uniform(dist) << "\n";
        //dist 1/3, 1/9, 1/9, 1/3, 1/9
        utils::print_distribution(dist);
    }

    void test_TVD()
    {
        vector<double> uniform;
        vector<double> dirac0;
        dirac0.push_back(1);
        for (int i = 0; i < 12; ++i)
        {
            uniform.push_back(1. / 12);
            if (i > 0)
                dirac0.push_back(0.);
        }
        cout << "DIRAC: ";
        for (auto i : dirac0)
            cout << i << ", ";
        cout << "\n";

        cout << "UNIFORM: ";
        for (auto i : uniform)
            cout << i << ", ";
        cout << "\n";

        cout << "Dirac 0: " << utils::get_TV_distance_from_uniform(dirac0) << "\n";
        cout << "Uniform: " << utils::get_TV_distance_from_uniform(uniform) << "\n";
    }

    void test_all_transpositions(int graph_size)
    {
        auto transps = utils::get_all_transpositions(graph_size);
        for (auto i : transps)
        {
            cout << i.first << " " << i.second << "\n";
        }
        cout << transps.size() << "\n";
    }

}

////////////////////////////////////////////////////

int main()
{
    cout << fixed << setprecision(10);

    int graph_size = 10000;

    simulation::properties CS;
    CS.steps = 100000;
    CS.RNG = utils::get_RNG();
    //CS.epsilon = 10.*numeric_limits<double>::epsilon();
    CS.epsilon = 0.0001;
    CS.fragmentation = false;

    CS.DEBUG = false;
    CS.coag_UOC = false;
    //CS.print_cycles = true;
    //CS.init_graph = init_graph::generate_random_permutation(graph_size, CS.RNG);
    CS.init_graph = init_graph::generate_singletons(graph_size);

    //simulation::get_averaged_TVD_quenched(CS, 500);
    simulation::simulate_quenched_biggest_ATM(CS);


}