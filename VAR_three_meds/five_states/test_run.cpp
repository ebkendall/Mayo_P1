#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

void generateStateSequences_start(const mat& adj_matrix, int current_state, int sequence_length, vec current_sequence, mat& result_sequences, int& sequence_index) {
    // Add current state to the sequence
    current_sequence(sequence_length - 1) = current_state;
    
    // Base case: If sequence_length is reached, add current_sequence to results
    if (sequence_length == 1) {
        result_sequences.row(sequence_index) = current_sequence.t();
        sequence_index++;
        return;
    }
    
    // Recursively explore all possible next states
    for (int next_state = 0; next_state < adj_matrix.n_rows; ++next_state) {
        if (adj_matrix(current_state, next_state) == 1) {
            generateStateSequences_start(adj_matrix, next_state, sequence_length - 1, current_sequence, result_sequences, sequence_index);
        }
    }
}

// [[Rcpp::export]]
arma::mat stateSequences_start(int starting_state, int sequence_length, const arma::mat& adjacency_matrix) {
    // Initialize result matrix to store all sequences
    int num_sequences = pow(adjacency_matrix.n_rows, sequence_length);
    mat result_sequences(num_sequences, sequence_length);
    
    // Initialize current sequence vector
    vec current_sequence(sequence_length);
    
    // Call recursive function to generate sequences
    int sequence_index = 0;
    generateStateSequences_start(adjacency_matrix, starting_state, sequence_length, current_sequence, result_sequences, sequence_index);
    
    // Resize result_sequences to actual number of sequences found
    result_sequences.resize(sequence_index, sequence_length);
    
    return result_sequences;
}


void generateStateSequences(const mat& adj_matrix, int sequence_length, vec current_sequence, mat& result_sequences, int end_state, int& sequence_index) {
    // Base case: If sequence_length is 1, check if we end at the end_state
    if (sequence_length == 1) {
        for (int i = 0; i < adj_matrix.n_rows; ++i) {
            if (adj_matrix(i, end_state) == 1) {
                // Add current sequence to results
                current_sequence(0) = i;
                result_sequences.row(sequence_index) = current_sequence.t();
                sequence_index++;
            }
        }
        return;
    }
    
    // Recursively explore all possible previous states
    for (int prev_state = 0; prev_state < adj_matrix.n_rows; ++prev_state) {
        if (adj_matrix(prev_state, end_state) == 1) {
            generateStateSequences(adj_matrix, sequence_length - 1, current_sequence, result_sequences, prev_state, sequence_index);
        }
    }
}

// [[Rcpp::export]]
arma::mat stateSequences(int end_state, int sequence_length, const arma::mat& adjacency_matrix) {
    // Initialize result matrix to store all sequences
    int num_sequences = pow(adjacency_matrix.n_rows, sequence_length);
    mat result_sequences(num_sequences, sequence_length);
    
    // Initialize current sequence vector
    vec current_sequence(sequence_length);
    
    // Call recursive function to generate sequences
    int sequence_index = 0;
    generateStateSequences(adjacency_matrix, sequence_length, current_sequence, result_sequences, end_state, sequence_index);
    
    // Resize result_sequences to actual number of sequences found
    result_sequences.resize(sequence_index, sequence_length);
    
    return result_sequences;
}

// // Recursive function to find all state sequences of given length ending in endState
// void findStateSequences_backward(int currentState,
//                         int endState,
//                         const arma::imat& adjacencyMatrix,
//                         int sequenceLength,
//                         std::vector<std::vector<int>>& allSequences,
//                         std::vector<int>& currentSequence) {
//     
//     // Add current state to the current sequence
//     currentSequence.push_back(currentState);
//     
//     // If current sequence length matches desired length and current state is endState, add to results
//     if (currentSequence.size() == sequenceLength && currentState == endState) {
//         allSequences.push_back(currentSequence);
//     }
//     
//     // If current sequence length is less than desired length, continue exploring
//     if (currentSequence.size() < sequenceLength) {
//         // Find next states
//         arma::uvec nextStates = find(adjacencyMatrix.row(currentState) == 1);
//         
//         // Recursively find sequences of desired length ending in endState starting from each next state
//         for (size_t i = 0; i < nextStates.n_elem; ++i) {
//             findStateSequences_backward(nextStates(i), endState, adjacencyMatrix, sequenceLength, allSequences, currentSequence);
//         }
//     }
//     
//     // Remove current state from the current sequence (backtracking)
//     currentSequence.pop_back();
// }
// 
// // [[Rcpp::export]]
// List getAllStateSequences_backward(int end_state, const arma::imat& adjacency_matrix, int sequence_length) {
//     int nStates = adjacency_matrix.n_rows;
//     
//     // Check validity of end_state
//     if (end_state < 0 || end_state >= nStates) {
//         stop("Invalid ending state index.");
//     }
//     
//     std::vector<std::vector<int>> allSequences;
//     std::vector<int> currentSequence;
//     
//     // Find all sequences of specified length ending in end_state
//     for (int start_state = 0; start_state < nStates; ++start_state) {
//         findStateSequences_backward(start_state, end_state, adjacency_matrix, sequence_length, allSequences, currentSequence);
//     }
//     
//     // // Convert C++ vectors to R list of integer vectors
//     List result(allSequences.size());
//     for (size_t i = 0; i < allSequences.size(); ++i) {
//         result[i] = wrap(allSequences[i]);
//     }
//     
//     return result;
// }
// 
// // Recursive function to find all state sequences of given length starting from startState
// void findStateSequences_forward(int currentState,
//                         const arma::imat& adjacencyMatrix,
//                         int sequenceLength,
//                         std::vector<std::vector<int>>& allSequences,
//                         std::vector<int>& currentSequence) {
//     
//     // Add current state to the current sequence
//     currentSequence.push_back(currentState);
//     
//     // If current sequence length matches desired length, add to results
//     if (currentSequence.size() == sequenceLength) {
//         allSequences.push_back(currentSequence);
//     }
//     
//     // If current sequence length is less than desired length, continue exploring
//     if (currentSequence.size() < sequenceLength) {
//         // Find next states
//         arma::uvec nextStates = find(adjacencyMatrix.row(currentState) == 1);
//         
//         // Recursively find sequences of desired length starting from each next state
//         for (size_t i = 0; i < nextStates.n_elem; ++i) {
//             findStateSequences_forward(nextStates(i), adjacencyMatrix, sequenceLength, allSequences, currentSequence);
//         }
//     }
//     
//     // Remove current state from the current sequence (backtracking)
//     currentSequence.pop_back();
// }
// 
// // [[Rcpp::export]]
// List getAllStateSequences_forward(int start_state, const arma::imat& adjacency_matrix, int sequence_length) {
//     int nStates = adjacency_matrix.n_rows;
//     
//     // Check validity of start_state
//     if (start_state < 0 || start_state >= nStates) {
//         stop("Invalid starting state index.");
//     }
//     
//     std::vector<std::vector<int>> allSequences;
//     std::vector<int> currentSequence;
//     
//     // Find all sequences of specified length starting from start_state
//     findStateSequences_forward(start_state, adjacency_matrix, sequence_length, allSequences, currentSequence);
//     
//     // Convert C++ vectors to R list of integer vectors
//     List result(allSequences.size());
//     for (size_t i = 0; i < allSequences.size(); ++i) {
//         result[i] = wrap(allSequences[i]);
//     }
//     
//     return result;
// }
// 
// void findStateSequences_both(int currentState,
//                         int endState,
//                         const arma::imat& adjacencyMatrix,
//                         int sequenceLength,
//                         std::vector<std::vector<int>>& allSequences,
//                         std::vector<int>& currentSequence) {
//     
//     // Add current state to the current sequence
//     currentSequence.push_back(currentState);
//     
//     // If current sequence length matches desired length and current state is endState, add to results
//     if (currentSequence.size() == sequenceLength && currentState == endState) {
//         allSequences.push_back(currentSequence);
//     }
//     
//     // If current sequence length is less than desired length, continue exploring
//     if (currentSequence.size() < sequenceLength) {
//         // Find next states
//         arma::uvec nextStates = find(adjacencyMatrix.row(currentState) == 1);
//         
//         // Recursively find sequences of desired length starting from each next state
//         for (size_t i = 0; i < nextStates.n_elem; ++i) {
//             findStateSequences_both(nextStates(i), endState, adjacencyMatrix, sequenceLength, allSequences, currentSequence);
//         }
//     }
//     
//     // Remove current state from the current sequence (backtracking)
//     currentSequence.pop_back();
// }
// 
// // [[Rcpp::export]]
// List getAllStateSequences_both(int start_state, int end_state, const arma::imat& adjacency_matrix, int sequence_length) {
//     int nStates = adjacency_matrix.n_rows;
//     
//     // Check validity of start_state and end_state
//     if (start_state < 0 || start_state >= nStates || end_state < 0 || end_state >= nStates) {
//         stop("Invalid starting or ending state index.");
//     }
//     
//     std::vector<std::vector<int>> allSequences;
//     std::vector<int> currentSequence;
//     
//     // Find all sequences of specified length starting from start_state and ending in end_state
//     findStateSequences_both(start_state, end_state, adjacency_matrix, sequence_length, allSequences, currentSequence);
//     
//     // Convert C++ vectors to R list of integer vectors
//     List result(allSequences.size());
//     for (size_t i = 0; i < allSequences.size(); ++i) {
//         result[i] = wrap(allSequences[i]);
//     }
//     
//     return result;
// }


// [[Rcpp::export]]
void test_fnc() {
    int nu_R = 1000;
    //   arma::mat Psi_R(4,4,arma::fill::eye);
    arma::vec scalar_vec_R = {4.58, 98.2, 101.3, 7.6};
    scalar_vec_R = (nu_R - 4 - 1) * scalar_vec_R;
    arma::mat Psi_R = arma::diagmat(scalar_vec_R);
    
    Rcpp::Rcout << Psi_R << std::endl;
}