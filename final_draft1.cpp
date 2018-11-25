#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <sys/time.h>
#include <string>
#include <ctime>
#include <map>
using namespace std;

std::vector<std::vector<int> > generate_random_test_sequence(int num_inputs, int num_vec){
    std::vector<std::vector<int> > test_sequence;
    std::vector<int> temp;
    for (int i =0; i < num_vec; i++){
        for (int j=0; j < num_inputs; j++){
            int a = rand()%2;
            temp.push_back(a);
        }
        test_sequence.push_back(temp);
        temp.clear();
    }
    return test_sequence;
}


void print_test_seq(std::vector<std::vector<int> > test_seq){
    for (int i =0; i < test_seq.size(); i++){
        for (int j=0; j< test_seq[i].size(); j++){
            cout<< test_seq[i][j] << " ";
        }
        cout<<"\n";
    }
}

// main()
main(int argc, char *argv[])
{
    int num_inputs = 12;//circuit -> numInputs;
    int num_ff = 12;//circuit -> numFF;
    int num_vectors = 10;
    std::vector<int> count_zero_ffs, count_one_ffs, ff_bias;

    for( int count = 0; count < num_ffs; count++){
        count_zero_ffs.push_back(0);
        count_one_ffs.push_back(0);
    }

///////////////////////////////////////////////////////////////////////////////////////////
//                    CONTROLLABILITY BASED PARTITIONING
///////////////////////////////////////////////////////////////////////////////////////////
    //Let number of inputs = n
    for(int count =0; count < 500; count++){
        std::vector<std::vector<int> > test_seq = generate_random_test_sequence( num_inputs, num_vectors);
        //print_test_seq(test_seq);
        //ff_states
        std::vector<int> ff_states = circuit->logicsim(test_seq);

        //count
        for(int i = 0; i < ff_states.size(); i++){
            if(ff_states[i] == 0){
                count_zero_ffs[i]++;
            }else{
                count_one_ffs[i]++;
            }
        }
    }

///////////////////////////////////////////////////////////////////////////////////////
//  ff_bias: ( 0 - 1 ) We will divide this into 5 state partitioning
//  state_partition vector stores the state number for every flip flop
//  0 - 0.2: S0     0.2 - 0.4: S1      0.4-0.6: S2      0.6-0.8: S3     0.8-1.0: S4
///////////////////////////////////////////////////////////////////////////////////////
    //Controllability Based Partitioning
    std::vector<int> state_partition;
    for(int i = 0; i < ff_states.size(); i++){
        ff_bias.push_back(abs(count_one_ffs[i] - count_zero_ffs[i])/500);

        if(ff_bias[i] <= 0.2){
            state_partition.push_back(0);
        }else if(ff_bias[i] <= 0.4){
            state_partition.push_back(1);
        }else if(ff_bias[i] <= 0.6){
            state_partition.push_back(2);
        }else if(ff_bias[i] <= 0.8){
            state_partition.push_back(3);
        }else if(ff_bias[i] <= 1){
            state_partition.push_back(4);
        }
    }
//////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< std::vector < std::vector < int> > > test_seq_population, curr_population, next_population;
    std::map< std::vector<int>, int > table_state0, table_state1, table_state2, table_state3, table_state4;


///////////////////////////////////////////////////////////////////////////////////////////////
//          DO GENETIC ALGORITHM 1000 TIMES
///////////////////////////////////////////////////////////////////////////////////////////////
    for( int times =0; times < 1000; times++){
        //initialize a random current population of 1000 vectors
        for( int a = 0; a<1000; a++) curr_population.push_back(generate_random_test_sequence(num_inputs, num_vectors));
        int num_test_seqs = curr_population.size();
        std::vector< int> fitness;
//////////////////////////////////////////////////////////////////////////////////////////////////////
//      COMPUTE FITNESS
//////////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector< int> curr_state_0, curr_state_1, curr_state_2, curr_state_3, curr_state_4;
        for( int i =0; i < num_test_seqs ; i++){
                std::vector<int> states = circuit->logicsim(curr_population[i]);      //Flip flop states
                int f = 0; //compute_fit( state_partition, states);          //compute fitness
                for(int j =0; j < states.size(); j++){
                    if(state_partition[j] == 0){
                        curr_state_0.push_back(states[j]);
                    }else if(state_partition[j] == 1){
                        curr_state_1.push_back(states[j]);
                    }else if(state_partition[j] == 2){
                        curr_state_2.push_back(states[j]);
                    }else if(state_partition[j] == 3{
                        curr_state_3.push_back(states[j]);
                    }else if(state_partition[j] == 4){
                        curr_state_4.push_back(states[j]);
                    }
                }

                std::map< std::vector<int>, int >::iterator it;

                if(table_state_0.count(curr_state_0)>0){
                    it = table_state_0.find(curr_state_0);
                    f += 0.2*(1/it->second);
                    it->second++;
                }else{
                    table_state_0[curr_state_0] = 1;
                }

                if(table_state_1.count(curr_state_1)>0){
                    it = table_state_1.find(curr_state_1);
                    f += 0.4*(1/it->second);
                    it->second++;
                }else{
                    table_state_1[curr_state_1] = 1;
                }

                if(table_state_2.count(curr_state_2)>0){
                    it = table_state_2.find(curr_state_2);
                    f += 0.6*(1/it->second);
                    it->second++;
                }else{
                    table_state_2[curr_state_2] = 1;
                }

                if(table_state_3.count(curr_state_3)>0){
                    it = table_state_3.find(curr_state_3);
                    f += 0.8*(1/it->second);
                    it->second++;
                }else{
                    table_state_3[curr_state_3] = 1;
                }

                if(table_state_4.count(curr_state_4)>0){
                    it = table_state_4.find(curr_state_4);
                    f += (1/it->second);
                    it->second++;
                }else{
                    table_state_4[curr_state_4] = 1;
                }
                curr_state_0.clear(); curr_state_1.clear(); curr_state_2.clear();
                curr_state_3.clear(); curr_state_4.clear();

                fitness.push_back(f);
        }
///////////////////////////////////////////////////////////////////////////////////////////////////////
//              FITNESS COMPUTED
///////////////////////////////////////////////////////////////////////////////////////////////////////
//      GENETIC ALGORITHM FOR max_gen GENERATIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////
        for( int gen_num = 0; gen_num < max_gen ; gen_num++){
            //Getting test sequence with best fitness
            int index = std::max_element(fitness.begin(),fitness.end()) - fitness.begin();
            test_seq_population.push_back(curr_population[index]);
            next_population.clear();

            for(int j = 1; j < curr_population.size()/2; j++){
                //////////////////////////////////////////////////////////////////////////////
                //              RANDOMLY SELECT TWO PARENTS
                //////////////////////////////////////////////////////////////////////////////
                random_par1 = (curr_population.size() + rand())%curr_population.size();
                random_par2 = (curr_population.size() + rand())%curr_population.size();

                std::vector< std::vector< int > > parent1 = curr_population[random_par1];
                std::vector< std::vector< int > > parent2 = curr_population[random_par2];

                std::vector< std::vector< int > > child1, child2;
                std::vector<int> temp1, temp2;
                ///////////////////////////////////////////////////////////////////////////////////////////
                //                     SELECTED TWO PARENTS
                ///////////////////////////////////////////////////////////////////////////////////////////
                //                      UNIFORM CROSSOVER
                ///////////////////////////////////////////////////////////////////////////////////////////
                // Uniform crossover( random_par1, random_par2, child1, child2)
                std::vector< std::vector< int > > random_mask = generate_random_test_sequence(num_inputs, num_vec);

                for(int vec =0; vec < num_vec ; vec++){
                    for(int in =0; in < num_inputs; in++){
                            if(random_mask[vec][in] == 0){
                                temp1.push_back(parent1[vec][in]);
                                temp2.push_back(parent2[vec][in]);
                            }else if(random_mask[vec][in] == 1){
                                temp1.push_back(parent2[vec][in]);
                                temp2.push_back(parent1[vec][in]);
                            }
                    }
                    child1.push_back(temp1); child2.push_back(temp2);
                    temp1.clear(); temp2.clear();
                }
                //////////////////////////////////////////////////////////////////////////////////////////
                //              CROSSOVER DONE
                ///////////////////////////////////////////////////////////////////////////////////////////
                //                  MUTATION
                ///////////////////////////////////////////////////////////////////////////////////////////
                //1% mutation of child1 and child2
                // To determine 1% mutation we will first calculate how many bits should we flip
                // Then we will randomly select said number of bits
                int num_bits_to_flip = 0.01*num_vec*num_inputs;
                if(num_bits_to_flip == 0) num_bits_to_flip++;

                int range = num_vec*num_inputs;
                for(int flips =0; flips < num_bits_to_flip; flips++){
                    int bit = rand()%range;
                    child1[bit/num_inputs][bit%num_inputs] = 1 - child1[bit/num_inputs][bit%num_inputs];
                    child2[bit/num_inputs][bit%num_inputs] = 1 - child2[bit/num_inputs][bit%num_inputs];
                }
                /////////////////////////////////////////////////////////////////////////////////////////
                //                  MUTATION DONE
                /////////////////////////////////////////////////////////////////////////////////////////

                next_population.push_back(child1); next_population.push_back(child2);
            }

            fitness.clear();
///////////////////////////////////////////////////////////////////////////////////////////
//          COMPUTE FITNESS
//////////////////////////////////////////////////////////////////////////////////////////
            for( int i =0; i < next_population.size() ; i++){
                std::vector<int> states = circuit->logicsim(next_population[i]);      //Flip flop states
                //int f = compute_fit( state_partition, states);          //compute fitness
                int f = 0; //compute_fit( state_partition, states);          //compute fitness
                for(int j =0; j < states.size(); j++){
                    if(state_partition[j] == 0){
                        curr_state_0.push_back(states[j]);
                    }else if(state_partition[j] == 1){
                        curr_state_1.push_back(states[j]);
                    }else if(state_partition[j] == 2){
                        curr_state_2.push_back(states[j]);
                    }else if(state_partition[j] == 3{
                        curr_state_3.push_back(states[j]);
                    }else if(state_partition[j] == 4){
                        curr_state_4.push_back(states[j]);
                    }
                }

                std::map< std::vector<int>, int >::iterator it;

                if(table_state_0.count(curr_state_0)>0){
                    it = table_state_0.find(curr_state_0);
                    f += 0.2*(1/it->second);
                    it->second++;
                }else{
                    table_state_0[curr_state_0] = 1;
                }

                if(table_state_1.count(curr_state_1)>0){
                    it = table_state_1.find(curr_state_1);
                    f += 0.4*(1/it->second);
                    it->second++;
                }else{
                    table_state_1[curr_state_1] = 1;
                }

                if(table_state_2.count(curr_state_2)>0){
                    it = table_state_2.find(curr_state_2);
                    f += 0.6*(1/it->second);
                    it->second++;
                }else{
                    table_state_2[curr_state_2] = 1;
                }

                if(table_state_3.count(curr_state_3)>0){
                    it = table_state_3.find(curr_state_3);
                    f += 0.8*(1/it->second);
                    it->second++;
                }else{
                    table_state_3[curr_state_3] = 1;
                }

                if(table_state_4.count(curr_state_4)>0){
                    it = table_state_4.find(curr_state_4);
                    f += (1/it->second);
                    it->second++;
                }else{
                    table_state_4[curr_state_4] = 1;
                }
                fitness.push_back(f);
            }
///////////////////////////////////////////////////////////////////////////////////////////
//          FITNESS COMPUTED
///////////////////////////////////////////////////////////////////////////////////////////

            curr_population = next_population;
            curr_state_0.clear(); curr_state_1.clear(); curr_state_2.clear();
            curr_state_3.clear(); curr_state_4.clear();
        }
////////////////////////////////////////////////////////////////////////////////////////////
//          GA DONE ONCE
////////////////////////////////////////////////////////////////////////////////////////////
    }
////////////////////////////////////////////////////////////////////////////////////////////
//          GA DONE 1000 TIMES
////////////////////////////////////////////////////////////////////////////////////////////
}


