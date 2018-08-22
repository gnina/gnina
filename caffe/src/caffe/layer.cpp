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
              LOG(ERROR) << "Backward relevance not implemented in layer " << this->type();
                abort(); // you need to implement
            }
        }
    }
}
INSTANTIATE_CLASS(Layer);

}  // namespace caffe
