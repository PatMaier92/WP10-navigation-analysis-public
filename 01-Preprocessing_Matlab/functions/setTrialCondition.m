function [TC]=setTrialCondition(condition,feedback)
% setTrialCondition Assigns correct trial condition for Starmaze WP10.
%
% Input: 
% condition is trial condition information (string). "main_ret" is missing,
% therefore this info needs to be corrected. 
% feedback indicates whether goal is visible (true/false) (string)
%
% Returns: TC is trial type (string)
% main_learn, main_ret, allo_ret, ego_ret, practise

training='main_learn';
% note: main_learn + no feedback (simple retrieval) is not marked as
% extra category in Starmaze WP10 output files yet.
allo='main_allo';
ego='main_ego';
mc='practise_motor'; 

TC_training=contains(condition,training);
TC_allo=contains(condition,allo);
TC_ego=contains(condition,ego);
TC_mc=contains(condition,mc); 

if TC_training==1
    if feedback=="false"
        TC="main_ret";
    else
        TC="main_learn";
    end
elseif TC_allo==1
    TC="allo_ret";
elseif TC_ego==1
    TC="ego_ret";
elseif TC_mc==1
    TC="practise";
else
    TC="999";
end

end
