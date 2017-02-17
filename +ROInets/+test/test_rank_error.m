classdef test_rank_error < matlab.unittest.TestCase

	% Check that rank deficient node timeseries are correctly identified
	% Depending on the orthogonalization method, rank errors can be thrown 
	% at different levels. They should be propagated through remove_source_leakage()

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

