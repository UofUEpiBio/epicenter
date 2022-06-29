#include "epiworld.hpp"
#include <iostream>

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
    ventilator,
    net_id
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

// Function to distribute the virus
std::vector<std::vector< Agent<int>* >> group_composition;

VirusToAgentFun<int> dist_virus_fun = [](Virus<int> & v, Model<int> * m) -> void {


    // First time it runs, we need to build the groups
    group_composition.clear();
    for (Agent<int> & a : (*m->get_agents()))
    {

        size_t id = static_cast<size_t>(a(CNAMES::net_id));
        if (group_composition.size() < (id + 1))
            group_composition.resize(id + 1);

        group_composition[id].push_back(&a);

    }
    
    /// Adding the virus, one per group
    for (auto & g : group_composition)
    {
        size_t id = static_cast<size_t>(std::floor(m->runif() * g.size()));
        g[id]->add_virus(v);
    }

};


int main(int argc, char* argv[]) 
{

    std::string ntype;
    size_t nnets;
    if (argc == 1)
    {

        ntype = "ergm";
        nnets = 1000;

    } else if (argc != 3)
        throw std::logic_error("You have to pass the number of networks to read in.");
    else {

        ntype = argv[2];
        nnets = std::strtoul(argv[1], nullptr, 0);

    }

    std::cout << "Simulation for " << ntype << ". Netcount: " << nnets << std::endl;

    // Function to save the run
    std::function<void(size_t,Model<>*)> saver = nullptr;
    if (ntype == "ergm")
        saver = make_save_run<>(
            "simoutbreak-data/%03lu-ergm-sim.csv",
            true, false, false, false, false, false, true, true
            );
    else if (ntype == "degseq")
        saver = make_save_run<>(
            "simoutbreak-data/%03lu-degseq-sim.csv",
            true, false, false, false, false, false, true, true
            );
    else
        throw std::logic_error("The requested type is not available.");

    // Individuals' covariates
    std::vector< double > data;

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


    // Making sure they run with different seeding
    std::vector< size_t > seeds(nnets, 0);

    for (size_t i = 1u; i <= nnets; ++i)
    {
        
        // Bones of the model
        Model<> model;

        if (i == 1)
        {
            model.seed(0);
            for (size_t j = 0; j < nnets; ++j)
                seeds[j] = static_cast<size_t>(
                    std::floor(
                        model.runif() * std::numeric_limits<size_t>::max()
                    ));
        }

        // Suppressing output
        model.verbose_off();

        // Setting the data
        model.set_agents_data(&data[0u], 4u);

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
        model.add_param(1.0/3.0, "rec rate hcw");
        model.add_param(1.0/7.0, "rec rate (vent x <=65)");
        model.add_param(1.0/14.0, "rec rate (vent x >65)");
        model.add_param(1.0/3.5, "rec rate (no vent x <=65)");
        model.add_param(1.0/7.0, "rec rate (no vent x >65)");

        model.add_param(.3, "prob infect");
        model.add_param(5, "incubation");

        // Processing the network
        char buff[100];

        if (ntype == "ergm")
            snprintf(buff, sizeof(buff), "networks/ergm-%04lu.txt", i);
        else
            snprintf(buff, sizeof(buff), "networks/degseq-%04lu.txt", i);
        
        std::string fn = buff;
        model.agents_from_adjlist(fn, 4191);

        Virus<> disease("A virus");
        disease.set_status(States::Exposed, States::Recovered, States::Deseased);
        disease.set_queue(QueueValues::OnlySelf, QueueValues::Everyone);

        disease.set_prob_death_fun(prob_death);
        disease.set_prob_recovery_fun(prob_rec);
        disease.set_prob_infecting(&model("prob infect"));

        model.add_virus_fun(disease, dist_virus_fun);

        // seeds are presetted during the first iteration
        model.init(100, seeds[i - 1]);

        model.run();

        saver(i, &model);

        // Only print the final one
        if (i == 10u)
            model.print();

    }

    

    return 0;

}