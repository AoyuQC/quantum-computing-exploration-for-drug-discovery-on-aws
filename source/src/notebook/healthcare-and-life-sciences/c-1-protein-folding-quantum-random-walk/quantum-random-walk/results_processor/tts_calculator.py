import math

class TTS_calculator:

    def __init__(self, default_value_tts, index_min_energy, step, precision_solution):

        self.default_value_tts = default_value_tts
        self.index_min_energy = index_min_energy
        self.step = step
        self.precision_solution = precision_solution

    def calculate_tts_from_probability_matrix(self, probabilities_matrix):

        p_t = 0
        # if the index of min energy calculated by psi 4 is in the results of metropolis, p_t is extracted
        # else, the p_t is set to a very small value close to 0 (not 0 to avoid inf values)
        if self.index_min_energy in probabilities_matrix.keys():
            p_t = probabilities_matrix[self.index_min_energy]
        else:
            p_t = 0

        result = 0
        # Result is the calculated TTS
        if p_t >= 1:
            result = 1
        elif p_t == 0:
            result = self.default_value_tts
        else:
            result = self.step * (math.log10(1-self.precision_solution)/(math.log10(1-p_t)))

        return result