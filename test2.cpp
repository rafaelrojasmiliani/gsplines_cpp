struct A{
    virtual int b() = 0;
};

struct Aauxiliar: public A{
    int b() override { return 0;}
};

struct B{
    double b();
};

struct aux: public Aauxiliar, public B{
using B::b;
};



int main()
{
    aux var;
    return 0;
}
