classdef test_rank_error < matlab.unittest.TestCase

	% Some unit tests to verify properties of the phase randomization process

	methods(Test)
		function test_rank_deficient(self)
			nodeData = [1 0 0 0; 1 0 0 0];
			
			self.verifyError(@() ROInets.remove_source_leakage(nodeData,'symmetric'),'ROInets:RankError')
			self.verifyError(@() ROInets.remove_source_leakage(nodeData,'closest'),'ROInets:RankError')
			self.verifyError(@() ROInets.remove_source_leakage(nodeData,'householder'),'ROInets:RankError')
		end	

		function test_rank_not_rank_deficient(self)
			nodeData = [1 0 0 0; 1 1 0 0];
			
			self.verifyWarningFree(@() ROInets.remove_source_leakage(nodeData,'symmetric'))
			self.verifyWarningFree(@() ROInets.remove_source_leakage(nodeData,'closest'))
			self.verifyWarningFree(@() ROInets.remove_source_leakage(nodeData,'householder'))
		end	

	end
end

