#include "caffe/layer.hpp"

namespace caffe {

template <typename Dtype>
void Layer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top,
                const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
    for (int i = 0; i < top.size(); ++i) {

        if(propagate_down[i]) {
            if(top[i]->count() == bottom[i]->count()) {
                //copy diff
                caffe_copy(top[i]->count(), top[i]->cpu_diff(), bottom[i]->mutable_cpu_diff());
            }
            else
            {
                std::cout << "aborting in default layer\n";
                abort(); // you need to implement
            }
        }
    }
}

template <typename Dtype>
void Layer<Dtype>::Backward_relevance_split(const vector<Blob<Dtype>*>& top,
                const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom,
                const int blob_to_propagate)
{
    Backward_relevance(top, propagate_down, bottom);
}

INSTANTIATE_CLASS(Layer);

}  // namespace caffe
