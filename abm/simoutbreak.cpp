#include "epiworld.hpp"

using namespace epiworld;

enum States {
    Susceptible,
    Exposed,
    Infected,
    Deseased,
    Recovered
};

enum CNAMES {
    is_hcw,
    age,
    ventilator
};

// Incubation time
EPI_NEW_UPDATEFUN(update_exposed, int)
{
    if (m->runif() < (1.0 / m->par("incubation")))
        p->change_status(States::Infected, QueueValues::Everyone);
}

// Probability of death
EPI_NEW_VIRUSFUN(prob_death, int)
{

    // Extracting variables
    if (p->operator()(CNAMES::is_hcw) > .5)
        return m->par("death rate hcw");

    // Ventilator is either 0 or 1, so .5 makes it good enough
    double pdeath;
    if (p->operator()(CNAMES::ventilator) > .5)
    {

        if (p->operator()(CNAMES::age) <= 65) 
            pdeath = m->par("death rate (vent x <=65)");
        else
            pdeath = m->par("death rate (vent x >65)");

    } else {

        if (p->operator()(CNAMES::age) <= 65) 
            pdeath = m->par("death rate (no vent x <=65)");
        else
            pdeath = m->par("death rate (no vent x >65)");
        
    }

    return pdeath;

}

EPI_NEW_VIRUSFUN(prob_rec, int)
{

    // Extracting variables
    if (p->operator()(CNAMES::is_hcw) > .5)
        return m->par("rec rate hcw");

    // Ventilator is either 0 or 1, so .5 makes it good enough
    double prec;
    if (p->operator()(CNAMES::ventilator) > .5)
    {

        if (p->operator()(CNAMES::age) <= 65) 
            prec = m->par("rec rate (vent x <=65)");
        else 
            prec = m->par("rec rate (vent x >65)");

    } else {

        if (p->operator()(CNAMES::age) <= 65) 
            prec = m->par("rec rate (no vent x <=65)");
        else 
            prec = m->par("rec rate (no vent x >65)");
        
    }

    return prec;

}

int main() 
{

    // Bones of the model
    Model<> model;

    // Individuals' covariates
    // size_t n = 1000u;
    std::vector< double > data;
    // double prop_hcw   = .2;
    // double prop_older = .7;
    // double prop_vent  = .1;

    // // HCW or not
    // for (size_t i = 0u; i < n; ++i)
    //     if (model.runif() < prop_hcw)
    //         data[i] = 1.0;

    // // Older or not
    // for (size_t i = 0u; i < n; ++i)
    //     if (model.runif() < prop_older)
    //         data[n + i] = 1.0;

    // // Ventilator
    // for (size_t i = 0u; i < n; ++i)
    //     if (model.runif() < prop_vent)
    //         data[2 * n + i] = 1.0;

    // Reading the individual level data
    std::ifstream file_x("actor_attributes.txt");
    double x;
    while (!file_x.eof())
    {
        // Capturing x
        file_x >> x;
        
        // Assigning the value
        data.push_back(x);
    }

    model.set_agents_data(&data[0u], 3u);

    // Model states(statuses)
    model.add_status("Susceptible", sampler::make_update_susceptible<>({States::Exposed}));
    model.add_status("Exposed", update_exposed);
    model.add_status("Infected", default_update_exposed<>);
    model.add_status("Deseased");
    model.add_status("Recovered");

    // Prob of death
    model.add_param(0.01, "death rate hcw");
    model.add_param(0.20, "death rate (vent x <=65)");
    model.add_param(0.40, "death rate (vent x >65)");
    model.add_param(0.02, "death rate (no vent x <=65)");
    model.add_param(0.10, "death rate (no vent x >65)");

    // Prob of recovery
    model.add_param(1.0/7.0, "rec rate hcw");
    model.add_param(1.0/14.0, "rec rate (vent x <=65)");
    model.add_param(1.0/21.0, "rec rate (vent x >65)");
    model.add_param(1.0/7.0, "rec rate (no vent x <=65)");
    model.add_param(1.0/14.0, "rec rate (no vent x >65)");

    model.add_param(.3, "prob infect");
    model.add_param(7, "incubation");

    // Processing the network
    // model.agents_smallworld(10000);
    model.agents_from_adjlist("networks/ergm-0001.txt", 4191);

    Virus<> disease("A virus");
    disease.set_status(States::Exposed, States::Recovered, States::Deseased);
    disease.set_queue(QueueValues::OnlySelf, QueueValues::Everyone);

    disease.set_prob_death_fun(prob_death);
    disease.set_prob_recovery_fun(prob_rec);
    disease.set_prob_infecting(&model("prob infect"));

    model.add_virus_n(disease, 10);

    model.init(100, 221);

    model.run();

    model.print();

    return 0;

}